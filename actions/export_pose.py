"""Export a docked Vina PDBQT pose to SDF for parametrization.

This action reads ``docking/result.json`` from the direct parent job and converts
the docked Vina PDBQT pose into an SDF file for later parametrization. Grubicy is
used for job and parent-file resolution, while JSON parsing remains local because
``docking/result.json`` and ``export/result.json`` are gbsa-grub stage manifests.
The chemistry conversion is delegated to ``gbsa_pipeline.docking`` so this
workflow layer only validates files and writes declared action products.
"""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path
from typing import Any

from gbsa_pipeline.docking import export_pdbqt_to_sdf, load_first_sdf_molecule
from grubicy import open_job_from_directory, parent_file

PARENT_MANIFEST = "docking/result.json"
EXPORTED_SDF_NAME = "dockligand_vina_out.sdf"
ERROR_LOG_NAME = "export_pose_error.log"
RUN_LOG_NAME = "export_pose.log"


def _read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON object from a manifest file."""
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
    """Write a deterministic JSON manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _write_text(path: Path, text: str) -> Path:
    """Write a text log and create the parent directory first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _required_path(data: dict[str, Any], key: str) -> Path:
    """Return a required manifest path as a resolved file path."""
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


def _optional_path_string(data: dict[str, Any], key: str) -> str | None:
    """Return an optional manifest path as a resolved string when available."""
    if key not in data:
        return None

    return str(_required_path(data, key))


def _read_parent_manifest(job: Any) -> dict[str, Any]:
    """Read the docking manifest from the direct Grubicy parent job."""
    manifest_path = parent_file(job, PARENT_MANIFEST, must_exist=True)
    return _read_json(manifest_path)


def main(directory: str | Path) -> None:
    """Export the docked PDBQT pose to SDF and write the export manifest."""
    job = open_job_from_directory(str(directory))
    parent = _read_parent_manifest(job)

    pose_pdbqt = _required_path(parent, "pose_pdbqt")
    ligand_sdf = _required_path(parent, "ligand_sdf")
    protein_pdb = _optional_path_string(parent, "protein_pdb")

    export_dir = Path(job.fn("export"))
    export_dir.mkdir(parents=True, exist_ok=True)

    run_log = export_dir / RUN_LOG_NAME
    output_sdf = export_dir / EXPORTED_SDF_NAME
    result_json = export_dir / "result.json"

    _write_text(
        run_log,
        "\n".join(
            [
                "export_pose action started",
                f"parent_manifest = {PARENT_MANIFEST}",
                f"pose_pdbqt = {pose_pdbqt}",
                f"ligand_sdf = {ligand_sdf}",
                f"output_sdf = {output_sdf}",
                "",
            ]
        ),
    )

    try:
        template_mol = load_first_sdf_molecule(ligand_sdf, remove_hs=False)
    except Exception as exc:
        raise RuntimeError(f"Failed to load ligand template SDF: {ligand_sdf}") from exc

    try:
        exported_sdf = export_pdbqt_to_sdf(
            pose_pdbqt,
            output_sdf,
            template_mol=template_mol,
            add_hydrogens_after_template=True,
        )
    except Exception as exc:
        raise RuntimeError(
            "Failed to export docked PDBQT pose to SDF.\n"
            f"PDBQT pose: {pose_pdbqt}\n"
            f"Ligand template SDF: {ligand_sdf}\n"
            f"Expected output SDF: {output_sdf}"
        ) from exc

    exported_sdf = Path(exported_sdf).resolve()

    if not exported_sdf.exists():
        raise RuntimeError(
            "PDBQT-to-SDF export returned without creating the expected SDF.\n"
            f"Returned path: {exported_sdf}\n"
            f"Expected output SDF: {output_sdf}"
        )

    manifest = {
        "docked_sdf": str(exported_sdf),
        "pose_pdbqt": str(pose_pdbqt),
        "ligand_sdf": str(ligand_sdf),
        "score": parent.get("score"),
        "rank": parent.get("rank"),
        "engine": parent.get("engine"),
        "export_dir": str(export_dir),
    }

    if protein_pdb is not None:
        manifest["protein_pdb"] = protein_pdb

    _write_json(result_json, manifest)

    _write_text(
        run_log,
        "\n".join(
            [
                "export_pose action finished successfully",
                f"pose_pdbqt = {pose_pdbqt}",
                f"ligand_sdf = {ligand_sdf}",
                f"protein_pdb = {protein_pdb}",
                f"docked_sdf = {exported_sdf}",
                f"result_json = {result_json}",
                "",
            ]
        ),
    )


def _cli(argv: list[str]) -> int:
    """Run the action and write a compact traceback log on failure."""
    if len(argv) != 2:
        print(
            "Usage error: expected exactly one argument: the grubicy job directory.",
            file=sys.stderr,
        )
        return 2

    directory = Path(argv[1])

    try:
        main(directory)
        return 0
    except Exception as exc:  # noqa: BLE001
        report = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        log_path = directory / "export" / ERROR_LOG_NAME

        try:
            _write_text(log_path, report)
        except Exception as log_exc:  # noqa: BLE001
            report += (
                "\nCould not write export_pose error log\n"
                f"Target log path: {log_path}\n"
                f"Logging exception: {type(log_exc).__name__}: {log_exc}\n"
            )

        print(report, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli(sys.argv))
