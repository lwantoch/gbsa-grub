"""Run AutoDock Vina for the prepared receptor and ligand.

This action is the docking step of the first gbsa-grub smoke workflow. It reads
``prepare/result.json`` from the direct parent job and uses the prepared ligand
and receptor PDBQT files written by ``prepare_inputs.py``. The actual docking
execution is delegated to ``gbsa_pipeline.docking.VinaEngine`` so the workflow
layer does not reimplement Vina command construction or subprocess logging. The
action writes stable declared products under ``docking/`` so downstream actions
can consume a docked PDBQT pose and a small JSON manifest.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Any

from gbsa_pipeline.docking import DockingBox, DockingRequest, VinaEngine
from grubicy import open_job_from_directory, parent_file

PARENT_MANIFEST = "prepare/result.json"
EXPECTED_POSE_NAME = "dockligand_vina_out.pdbqt"
EXPECTED_LOG_NAME = "dockligand_vina.log"


def _read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON manifest from disk.

    The docking action consumes a lightweight manifest from the prepare step
    rather than guessing upstream paths. The decoded value must be a JSON object
    because downstream code expects named fields such as ``ligand_pdbqt``,
    ``receptor_pdbqt``, ``box_center``, and ``box_size``. This helper performs
    only generic JSON validation; action-specific validation happens in
    ``main``. The manifest intentionally stores file paths and metadata, not
    molecular structures.
    """
    manifest_path = Path(path)

    if not manifest_path.exists():
        raise FileNotFoundError(f"JSON manifest not found: {manifest_path}")

    if not manifest_path.is_file():
        raise ValueError(f"JSON manifest path is not a file: {manifest_path}")

    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in manifest: {manifest_path}") from exc

    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object in manifest: {manifest_path}")

    return data


def _write_json(path: Path, data: dict[str, Any]) -> Path:
    """Write a small JSON manifest with deterministic formatting.

    The docking action writes ``docking/result.json`` as the contract consumed
    by the pose-export action. The manifest records the final normalized pose
    path, the docking log path, the Vina score, and the prepared input files
    used for the run. Serialization errors are allowed to surface because a
    non-JSON-compatible manifest value is a runner bug. Scientific artefacts are
    written as separate declared products, not embedded into the JSON file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _read_parent_manifest(job: Any) -> dict[str, Any]:
    """Read ``prepare/result.json`` from the direct parent job.

    Grubicy owns parent-job resolution through the dependency declared in
    ``pipeline.toml``. This helper therefore asks Grubicy for the parent file
    instead of reconstructing parent paths manually from signac details. JSON
    parsing remains local because this file is a gbsa-grub stage manifest, not a
    Grubicy action-parameter model. Missing parent products fail before Vina is
    launched.
    """
    manifest_path = parent_file(job, PARENT_MANIFEST, must_exist=True)
    return _read_json(manifest_path)


def _required_path(data: dict[str, Any], key: str) -> Path:
    """Return a required manifest path as a resolved file path.

    The docking action needs concrete prepared input files from the parent
    prepare step. This helper keeps required-key, value-type, existence, and
    file checks explicit at the workflow boundary. Invalid manifests fail before
    the library docking layer is called, which gives a clearer runner-level
    error. The function does not inspect molecular contents.
    """
    if key not in data:
        available = ", ".join(sorted(data)) or "<none>"
        raise KeyError(
            f"Parent manifest is missing required key: {key}\n"
            f"Required parent manifest: {PARENT_MANIFEST}\n"
            f"Available keys: {available}"
        )

    value = data[key]
    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            f"Parent manifest key '{key}' must be a non-empty string path. "
            f"Got: {value!r}"
        )

    path = Path(value).resolve()

    if not path.exists():
        raise FileNotFoundError(f"Required parent file not found for '{key}': {path}")

    if not path.is_file():
        raise ValueError(f"Required parent path for '{key}' is not a file: {path}")

    return path


def _triple_from_manifest(data: dict[str, Any], key: str) -> tuple[float, float, float]:
    """Convert a manifest box field to a three-float tuple.

    The prepare action stores docking box values as JSON arrays. The library
    ``DockingBox`` expects fixed-length tuples of floats, so the conversion is
    explicit at the workflow boundary. Invalid values fail before Vina is
    launched, which makes malformed manifests easier to diagnose. The helper
    does not interpret the scientific meaning of the box; it only validates the
    expected shape of the data.
    """
    if key not in data:
        available = ", ".join(sorted(data)) or "<none>"
        raise KeyError(
            f"Parent manifest is missing required key: {key}\n"
            f"Required parent manifest: {PARENT_MANIFEST}\n"
            f"Available keys: {available}"
        )

    value = data[key]
    if not isinstance(value, list | tuple):
        raise ValueError(
            f"Parent manifest key '{key}' must be a sequence of three numbers. "
            f"Got: {value!r}"
        )

    if len(value) != 3:
        raise ValueError(
            f"Parent manifest key '{key}' must contain exactly three numbers. "
            f"Got: {value!r}"
        )

    try:
        return tuple(float(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Parent manifest key '{key}' must contain only numeric values. "
            f"Got: {value!r}"
        ) from exc


def _move_product_if_needed(source: Path, target: Path) -> Path:
    """Move a generated product to the declared workflow product path.

    ``VinaEngine`` names output files from the input ligand stem. The grubicy
    workflow declares stable ``dockligand_*`` products so later actions do not
    depend on the exact ligand-preparation filename. If the generated path is
    already the declared path, no filesystem operation is needed. Otherwise the
    generated product is moved to the declared path, replacing a stale target
    from an earlier failed run if necessary. This keeps the action output names
    stable without changing the library-level Vina wrapper.
    """
    source = source.resolve()
    target = target.resolve()

    if source == target:
        return target

    if not source.exists():
        raise FileNotFoundError(f"Generated docking product not found: {source}")

    target.parent.mkdir(parents=True, exist_ok=True)

    if target.exists():
        target.unlink()

    source.replace(target)
    return target


def main(directory: str | Path) -> None:
    """Run Vina docking and write the declared docking products.

    The action opens the current job from the workspace directory supplied by
    row/grubicy. It resolves the direct parent manifest through Grubicy rather
    than through direct signac path reconstruction. Prepared input paths and box
    values are validated at this workflow boundary before ``VinaEngine`` is
    called. The function writes stable declared products under ``docking/`` and
    records a compact manifest for the downstream pose-export action.
    """
    job = open_job_from_directory(str(directory))

    if shutil.which("vina") is None:
        raise RuntimeError("vina not available in PATH.")

    parent = _read_parent_manifest(job)

    ligand_pdbqt = _required_path(parent, "ligand_pdbqt")
    receptor_pdbqt = _required_path(parent, "receptor_pdbqt")
    protein_pdb = _required_path(parent, "protein_pdb")
    ligand_sdf = _required_path(parent, "ligand_sdf")

    docking_dir = Path(job.fn("docking"))
    docking_dir.mkdir(parents=True, exist_ok=True)

    expected_pose = docking_dir / EXPECTED_POSE_NAME
    expected_log = docking_dir / EXPECTED_LOG_NAME
    result_json = docking_dir / "result.json"

    box = DockingBox(
        center=_triple_from_manifest(parent, "box_center"),
        size=_triple_from_manifest(parent, "box_size"),
    )

    request = DockingRequest(
        receptor=receptor_pdbqt,
        ligands=[ligand_pdbqt],
        box=box,
        workdir=docking_dir,
    )

    docking_result = VinaEngine(binary="vina").dock(request=request)

    if len(docking_result.poses) != 1:
        raise RuntimeError(
            f"Expected exactly one docked pose, got {len(docking_result.poses)}."
        )

    pose = docking_result.poses[0]

    if pose.metadata.get("returncode") != 0:
        raise RuntimeError(
            "Vina docking failed. "
            f"Return code: {pose.metadata.get('returncode')}. "
            f"Log: {pose.metadata.get('log_file')}"
        )

    generated_pose = Path(pose.pose_path)
    generated_log = Path(
        pose.metadata.get(
            "log_file",
            docking_dir / f"{ligand_pdbqt.stem}_vina.log",
        )
    )

    final_pose = _move_product_if_needed(generated_pose, expected_pose)
    final_log = _move_product_if_needed(generated_log, expected_log)

    _write_json(
        result_json,
        {
            "engine": docking_result.engine,
            "score": pose.score,
            "rank": pose.rank,
            "pose_pdbqt": str(final_pose),
            "vina_log": str(final_log),
            "protein_pdb": str(protein_pdb),
            "ligand_sdf": str(ligand_sdf),
            "ligand_pdbqt": str(ligand_pdbqt),
            "receptor_pdbqt": str(receptor_pdbqt),
            "box_center": list(box.center),
            "box_size": list(box.size),
            "docking_dir": str(docking_dir),
        },
    )


if __name__ == "__main__":
    main(sys.argv[1])
