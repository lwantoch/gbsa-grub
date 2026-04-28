"""Shared helpers for grubicy action runners.

This module is intentionally kept under ``actions`` because it is workflow glue,
not gbsa-pipeline library code. The helpers make all action runners use the same
job opening, parent-manifest reading, JSON writing, and BioSimSpace disk handoff
logic. Every MD action receives a parent ``system.gro`` and ``system.top`` from
``result.json``, loads those files into BioSimSpace, runs one stage, and writes
a new ``system.gro``/``system.top`` pair for the next action. This avoids passing
Python objects between grubicy jobs, which would not work reliably across row
submissions.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import signac

REPO_ROOT = Path.cwd()


def open_job_from_directory(directory: str | Path) -> Any:
    """Open the signac job represented by a grubicy/row workspace directory.

    Grubicy action runners receive the workspace directory as the first command
    line argument. In a normal signac project the workspace directory name is the
    job id, so this helper opens the existing job instead of creating a new one.
    The return type is intentionally ``Any`` because this script should stay a
    thin workflow adapter and should not introduce local signac protocols. If
    the directory is not a valid signac job workspace, signac will raise a clear
    error before any stage logic is executed.
    """
    project = signac.get_project()
    return project.open_job(id=Path(directory).name)


def job_dir(job: Any) -> Path:
    """Return the absolute workspace path for a signac job.

    All stage output directories are created below this path. The helper keeps
    the runners independent from signac's internal string-returning APIs and
    makes path composition explicit. It does not create the job directory,
    because signac/grubicy should already have created the workspace before the
    runner starts.
    """
    return Path(job.workspace()).resolve()


def stage_dir(job: Any, stage: str) -> Path:
    """Create and return the output directory for one action stage.

    Each action writes into a directory named after the stage, for example
    ``minimization_sd`` or ``equilibration_npt``. Existing directories are
    accepted so failed jobs can be re-run without deleting the whole signac
    workspace. Individual output files are overwritten by the functions that
    create them.
    """
    path = job_dir(job) / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


def statepoint(job: Any) -> dict[str, Any]:
    """Return the job state point as a plain dictionary.

    Signac exposes state points as mapping-like objects. A plain dictionary is
    easier to pass into Pydantic models, JSON manifests, and small local helper
    functions. The returned dictionary is a copy and mutating it will not update
    the signac job.
    """
    return dict(job.sp)


def action_name(job: Any, fallback: str | None = None) -> str:
    """Return the grubicy action name for a job.

    Grubicy normally stores the action name in the state point as ``action``.
    During manual testing it is useful to pass a hard-coded fallback from the
    runner filename. The function raises if neither is available because the
    action name also determines the stage output directory and manifest path.
    """
    sp = statepoint(job)
    value = sp.get("action") or fallback
    if not value:
        raise KeyError("Could not determine action name from job.sp['action'].")
    return str(value)


def parent_action_name(job: Any, parent: Any) -> str:
    """Return the stage directory name used by the parent job.

    Child jobs created by grubicy commonly carry ``parent_action`` in their own
    state point. If that key is not present, this helper falls back to the
    parent job's ``action`` key. This makes the runners tolerant to small
    differences in grubicy specs while still enforcing that every parent output
    is stored under ``<parent_action>/result.json``.
    """
    child_sp = statepoint(job)
    parent_sp = statepoint(parent)

    value = child_sp.get("parent_action") or parent_sp.get("action")
    if not value:
        raise KeyError("Could not determine parent action name.")

    return str(value)


def resolve_repo_path(value: str | Path) -> Path:
    """Resolve an input path from the state point.

    Absolute paths are used unchanged. Relative paths are interpreted relative to
    the repository root, which is expected to be the working directory when row
    executes the action command. The function does not check existence by
    itself because the caller can provide a more specific error depending on
    whether the path is a protein, ligand, index file, or optional input.
    """
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (REPO_ROOT / path).resolve()


def require_sp_path(job: Any, key: str) -> Path:
    """Read a required state point path and resolve it.

    The path key is required for source inputs such as ``protein_pdb`` and
    ``ligand_sdf``. A clear ``KeyError`` is more useful than a later failure in
    RDKit, Meeko, OpenMM, or BioSimSpace. The function resolves relative paths
    against the repository root.
    """
    sp = statepoint(job)
    if key not in sp:
        raise KeyError(f"Missing required state point key: {key}")
    return resolve_repo_path(sp[key])


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object from disk.

    Stage handoff uses small JSON manifests, always as dictionaries. Returning a
    dictionary avoids repeating type checks in every runner. If the file contains
    anything other than a JSON object, this helper raises a targeted error.
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise TypeError(f"Expected JSON object in {path}, got {type(data).__name__}.")
    return data


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Write a deterministic JSON manifest.

    The manifest stores paths, stage names, parameters, and selected native
    output files. Paths and Pydantic-like objects are converted to serializable
    values through ``_json_default``. Sorting keys and using indentation keeps
    the manifests inspectable during early grubicy debugging.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if hasattr(value, "value"):
        return value.value
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable.")


def parent_manifest(job: Any) -> tuple[Any, Path, dict[str, Any]]:
    """Return the parent job, parent manifest path, and decoded manifest.

    The child action reads its input system from the parent stage manifest.
    The manifest convention is always ``<parent_action>/result.json`` inside
    the parent workspace. The helper uses grubicy's parent accessors and fails
    before simulation starts if the declared parent output is missing.
    """
    from grubicy import get_parent, parent_path, parent_product_exists  # noqa: PLC0415

    parent = get_parent(job)
    parent_stage = parent_action_name(job, parent)
    relpath = Path(parent_stage) / "result.json"

    if not parent_product_exists(job, relpath):
        raise FileNotFoundError(f"Missing parent product: {parent_path(job) / relpath}")

    result_path = parent_path(job) / relpath
    return parent, result_path, read_json(result_path)


def parent_system_paths(job: Any) -> tuple[Path, Path, Path, dict[str, Any]]:
    """Return parent ``system.gro``/``system.top`` paths plus manifest metadata.

    MD actions only need two files from the parent stage: coordinates and
    topology. The setup and MD runners write both ``output_*`` and ``system_*``
    keys, so this helper accepts either naming convention. It validates that the
    files exist before BioSimSpace is called, which makes grubicy failures easier
    to diagnose.
    """
    _parent, result_path, manifest = parent_manifest(job)

    gro_value = manifest.get("output_gro") or manifest.get("system_gro")
    top_value = manifest.get("output_top") or manifest.get("system_top")

    if gro_value is None or top_value is None:
        raise KeyError(
            "Parent manifest must contain output_gro/output_top or system_gro/system_top."
        )

    gro = Path(str(gro_value))
    top = Path(str(top_value))

    if not gro.exists():
        raise FileNotFoundError(f"Parent GRO file not found: {gro}")
    if not top.exists():
        raise FileNotFoundError(f"Parent TOP file not found: {top}")

    return gro, top, result_path, manifest


def load_bss_system(gro_file: Path, top_file: Path) -> Any:
    """Load a GROMACS coordinate/topology pair as a BioSimSpace system.

    This is the disk-based replacement for passing a Python system object from
    one MD stage to the next. It uses BioSimSpace's normal GROMACS file loading
    path and performs no chemistry changes. All scientific validation remains in
    the existing gbsa-pipeline helpers and in GROMACS itself.
    """
    import BioSimSpace as BSS  # noqa: PLC0415

    return BSS.IO.readMolecules([str(gro_file), str(top_file)])


def save_bss_system(system: Any, output_base: Path) -> tuple[Path, Path]:
    """Save a BioSimSpace system as ``output_base.gro`` and ``output_base.top``.

    The output files are deleted first so stale files cannot be mistaken for a
    successful write. BioSimSpace is then asked to write GROMACS coordinate and
    topology formats. The expected files are checked explicitly because the next
    grubicy action depends on them.
    """
    import BioSimSpace as BSS  # noqa: PLC0415

    output_base.parent.mkdir(parents=True, exist_ok=True)

    output_gro = output_base.with_suffix(".gro")
    output_top = output_base.with_suffix(".top")

    for path in (output_gro, output_top):
        if path.exists():
            path.unlink()

    BSS.IO.saveMolecules(str(output_base), system, ["gro87", "grotop"])

    if not output_gro.exists():
        raise RuntimeError(f"Expected GRO output was not written: {output_gro}")
    if not output_top.exists():
        raise RuntimeError(f"Expected TOP output was not written: {output_top}")

    return output_gro, output_top


def collect_native_outputs(path: Path) -> dict[str, str]:
    """Collect useful native GROMACS files produced by one stage.

    BioSimSpace/GROMACS working directories may contain files such as ``.tpr``,
    ``.xtc``, ``.trr``, ``.edr``, ``.cpt``, ``.log``, ``.mdp``, and ``.ndx``.
    This helper stores the newest file found for each suffix. It does not rename
    or delete anything; it only exposes paths in ``result.json`` for later
    debugging and, eventually, GBSA.
    """
    suffixes = {".tpr", ".xtc", ".trr", ".edr", ".cpt", ".log", ".mdp", ".ndx"}
    found: dict[str, Path] = {}

    for candidate in path.rglob("*"):
        if not candidate.is_file():
            continue

        suffix = candidate.suffix.lower()
        if suffix not in suffixes:
            continue

        previous = found.get(suffix)
        if previous is None or candidate.stat().st_mtime > previous.stat().st_mtime:
            found[suffix] = candidate.resolve()

    return {
        suffix.lstrip("."): str(file_path)
        for suffix, file_path in sorted(found.items())
    }


def write_stage_manifest(
    *,
    job: Any,
    stage: str,
    parent_result: Path,
    input_gro: Path,
    input_top: Path,
    output_gro: Path,
    output_top: Path,
    params: dict[str, Any],
) -> None:
    """Write the standard result manifest for one MD action.

    Every MD action uses the same manifest keys so the next action can load its
    parent without knowing which exact stage produced the files. ``native_outputs``
    is included for debugging and for later GBSA integration. The manifest is
    written to ``<stage>/result.json`` in the current job workspace.
    """
    current_dir = stage_dir(job, stage)
    native_outputs = collect_native_outputs(current_dir)

    write_json(
        current_dir / "result.json",
        {
            "stage": stage,
            "parent_result": parent_result,
            "input_gro": input_gro,
            "input_top": input_top,
            "output_gro": output_gro,
            "output_top": output_top,
            "system_gro": output_gro,
            "system_top": output_top,
            "work_dir": current_dir,
            "params": params,
            "native_outputs": native_outputs,
            "tpr": native_outputs.get("tpr"),
            "trajectory_xtc": native_outputs.get("xtc"),
            "trajectory_trr": native_outputs.get("trr"),
            "edr": native_outputs.get("edr"),
            "log": native_outputs.get("log"),
            "index": native_outputs.get("ndx"),
        },
    )
