# actions/export_pose.py

"""Export a docked Vina PDBQT pose to SDF for parametrization.

This action is the pose-export step of the first gbsa-grub smoke workflow. It
reads ``docking/result.json`` from the direct parent job and converts the docked
PDBQT pose into an SDF file suitable for the later setup step. The actual ligand
chemistry reconstruction and PDBQT-to-SDF conversion are delegated to
``gbsa_pipeline.docking.export_pdbqt_to_sdf`` so this workflow layer does not
duplicate chemistry logic from the library. The action writes the declared
``export/dockligand_vina_out.sdf`` product plus a small ``export/result.json``
manifest for downstream actions, and writes a diagnostic error log when the
runner fails under row/grubicy.
"""

from __future__ import annotations

import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import signac
from gbsa_pipeline.docking import export_pdbqt_to_sdf, load_first_sdf_molecule
from grubicy import parent_path, parent_product_exists

PARENT_MANIFEST = "docking/result.json"
EXPORTED_SDF_NAME = "dockligand_vina_out.sdf"
ERROR_LOG_NAME = "export_pose_error.log"
RUN_LOG_NAME = "export_pose.log"


def _open_job(directory: str | Path) -> Any:
    """Open the signac job represented by the grubicy workspace directory.

    Grubicy passes the current job workspace directory to the action runner as
    the first positional argument. The directory basename is the signac job id,
    so the current project can reopen the job directly without a separate
    statepoint lookup. This keeps the action aligned with the grubicy runner
    model while avoiding a shared helper module for now. The returned job is
    used for parent-product checks and output path construction through
    ``job.fn(...)``.
    """
    return signac.get_project().open_job(id=Path(directory).name)


def _read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON manifest from disk.

    The export action consumes the docking-stage manifest rather than
    reconstructing upstream paths manually. The decoded JSON value must be an
    object because the runner expects named fields such as ``pose_pdbqt`` and
    ``ligand_sdf``. This helper performs only generic JSON validation while the
    required action-specific keys are checked in ``main``. Molecular data remain
    in normal workflow files and are not embedded in the manifest.
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

    The export action writes ``export/result.json`` as the contract consumed by
    the system-setup action. The manifest records the exported docked SDF, the
    source docked PDBQT pose, the original ligand SDF template, and available
    docking metadata. Serialization errors are allowed to surface because a
    non-JSON-compatible manifest value is a workflow bug. Large molecular
    artefacts are written as files, not embedded in this JSON.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _write_text(path: Path, text: str) -> Path:
    """Write diagnostic text and create the parent directory first.

    Row can hide the useful Python traceback when an action fails, so this
    action writes its own small text logs inside the job workspace. The helper is
    intentionally local to this runner and does not configure global logging.
    It is used for both the normal run log and the failure report. Filesystem
    errors are not suppressed during normal calls because missing logs should be
    visible while debugging the workflow.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _safe_json(data: Any) -> str:
    """Return a readable JSON representation for diagnostics.

    Diagnostic reports should be readable even when signac statepoints contain
    objects that are not directly JSON serializable. This helper uses
    ``default=str`` so the error report can still be written. It is used only
    for human-facing debug context and not for workflow manifests. Stable
    machine-readable manifests are still written through ``_write_json``.
    """
    return json.dumps(data, indent=2, sort_keys=True, default=str)


def _safe_directory_listing(path: Path, *, max_entries: int = 80) -> str:
    """Return a compact recursive listing for diagnostics.

    Failed grubicy actions are usually debugged by inspecting which declared
    products and parent products exist. This helper lists a bounded number of
    files below a directory so the error log shows whether the expected parent
    manifest, docked pose, or exported SDF was present. Missing directories are
    reported explicitly instead of raising. The function is diagnostic only and
    does not influence workflow behavior.
    """
    if not path.exists():
        return f"<missing directory: {path}>"

    if not path.is_dir():
        return f"<not a directory: {path}>"

    entries = sorted(item for item in path.rglob("*") if item.is_file())
    if not entries:
        return "<no files>"

    shown = entries[:max_entries]
    lines = [str(item.relative_to(path)) for item in shown]

    if len(entries) > max_entries:
        lines.append(f"... {len(entries) - max_entries} more files")

    return "\n".join(lines)


def _error_log_path(directory: str | Path | None) -> Path:
    """Return the best available error-log path for this action.

    Most failures happen after grubicy has created the action workspace, so the
    preferred location is ``<job>/export/export_pose_error.log``. If the CLI
    argument is missing or malformed, the fallback is a local
    ``export_pose_error.log`` in the current working directory. This makes the
    failure visible even when the runner cannot open the signac job. The helper
    does not require signac and therefore works for early CLI errors.
    """
    if directory is None:
        return Path(ERROR_LOG_NAME).resolve()

    return Path(directory).resolve() / "export" / ERROR_LOG_NAME


def _build_failure_report(
    *,
    directory: str | Path | None,
    exc: BaseException,
) -> str:
    """Build a detailed failure report for row/grubicy action errors.

    Row often reports only that an action exited with a non-zero code. This
    report records the exception type, full traceback, action workspace, current
    statepoint if readable, direct parent path if resolvable, parent manifest
    content if available, and a compact file listing. The report is intentionally
    broad because failures can come from missing parent products, malformed JSON,
    missing manifest keys, missing molecular files, RDKit/Meeko conversion
    errors, or output-write errors. The function must not raise while building
    diagnostics, so each optional section is protected independently.
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    directory_path = Path(directory).resolve() if directory is not None else None

    parts = [
        "export_pose action failed",
        "=========================",
        "",
        f"UTC time: {timestamp}",
        f"Exception type: {type(exc).__name__}",
        f"Exception message: {exc}",
        "",
        "Traceback:",
        "----------",
        "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)),
    ]

    if directory_path is None:
        parts.extend(["", "Action directory: <missing CLI argument>"])
        return "\n".join(parts)

    parts.extend(["", f"Action directory: {directory_path}"])

    statepoint_path = directory_path / "signac_statepoint.json"
    try:
        parts.extend(
            [
                "",
                f"Statepoint path: {statepoint_path}",
                "Statepoint:",
                "-----------",
                _safe_json(_read_json(statepoint_path)),
            ]
        )
    except Exception as statepoint_exc:  # noqa: BLE001
        parts.extend(
            [
                "",
                f"Statepoint path: {statepoint_path}",
                f"Could not read statepoint: {type(statepoint_exc).__name__}: {statepoint_exc}",
            ]
        )

    try:
        job = _open_job(directory_path)
        parent_dir = parent_path(job)
        parent_manifest_path = parent_dir / PARENT_MANIFEST

        parts.extend(
            [
                "",
                f"Resolved parent directory: {parent_dir}",
                f"Expected parent manifest: {parent_manifest_path}",
                f"Parent product exists according to grubicy: {parent_product_exists(job, PARENT_MANIFEST)}",
            ]
        )

        if parent_manifest_path.exists():
            parts.extend(
                [
                    "",
                    "Parent manifest:",
                    "----------------",
                    _safe_json(_read_json(parent_manifest_path)),
                ]
            )
        else:
            parts.extend(["", "Parent manifest: <missing>"])

        parts.extend(
            [
                "",
                "Parent directory file listing:",
                "------------------------------",
                _safe_directory_listing(parent_dir),
            ]
        )
    except Exception as parent_exc:  # noqa: BLE001
        parts.extend(
            [
                "",
                f"Could not resolve/read parent context: {type(parent_exc).__name__}: {parent_exc}",
            ]
        )

    parts.extend(
        [
            "",
            "Action directory file listing:",
            "------------------------------",
            _safe_directory_listing(directory_path),
            "",
            "Required parent manifest keys for this action:",
            "----------------------------------------------",
            "- pose_pdbqt",
            "- ligand_sdf",
            "",
            "Common fix if ligand_sdf is missing:",
            "-----------------------------------",
            'actions/dock.py must propagate "ligand_sdf" from prepare/result.json into docking/result.json.',
            'setup_system.py will also need "protein_pdb" later, so propagate that too.',
        ]
    )

    return "\n".join(parts)


def _read_parent_manifest(job: Any) -> dict[str, Any]:
    """Read ``docking/result.json`` from the direct parent job.

    Grubicy tracks the parent job through the dependency declared in
    ``pipeline.toml``. This helper checks the declared parent product before
    reading it, so the action follows the same completion model used by
    ``grubicy status`` and ``grubicy submit``. Missing parent products fail
    before any ligand conversion starts. The returned dictionary is interpreted
    by ``main`` because the required fields are action-specific.
    """
    if not parent_product_exists(job, PARENT_MANIFEST):
        parent_dir = parent_path(job)
        raise FileNotFoundError(
            "Missing parent product required by export_pose.\n"
            f"Required product: {PARENT_MANIFEST}\n"
            f"Resolved parent directory: {parent_dir}\n"
            f"Expected path: {parent_dir / PARENT_MANIFEST}"
        )

    return _read_json(parent_path(job) / PARENT_MANIFEST)


def _required_path(data: dict[str, Any], key: str) -> Path:
    """Return a required manifest path as a resolved ``Path``.

    The export action needs concrete files from the docking and prepare stages:
    the docked PDBQT pose and the original ligand SDF template. This helper keeps
    required-key and file-existence checks explicit at the workflow boundary.
    Missing keys and missing files should fail before the library conversion is
    called, giving a clearer action-level error. The function does not parse or
    validate molecular contents.
    """
    if key not in data:
        available = ", ".join(sorted(data)) or "<none>"
        raise KeyError(
            f"Parent manifest is missing required key: {key}\n"
            f"Required parent manifest: {PARENT_MANIFEST}\n"
            f"Available keys: {available}\n"
            "For this workflow, docking/result.json must contain at least "
            '"pose_pdbqt" and "ligand_sdf".'
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
    """Return an optional manifest path as a resolved string when available.

    The export result should propagate ``protein_pdb`` when the docking manifest
    provides it, because the setup action needs the original protein file. Older
    docking manifests may not contain that key while the workflow is being
    debugged, so this helper keeps optional propagation separate from required
    export inputs. If the key is present, the path is checked so stale or broken
    manifests are still reported clearly. ``None`` means the key was absent.
    """
    if key not in data:
        return None

    path = _required_path(data, key)
    return str(path)


def main(directory: str | Path) -> None:
    """Export the docked PDBQT pose to SDF and write the export manifest."""
    job = _open_job(directory)
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
                f"parent_manifest = {PARENT_MANIFEST}",
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
    """CLI wrapper that writes a failure report before returning a non-zero code.

    Unit tests should call ``main(...)`` directly so exceptions remain normal
    pytest failures. The command-line entry point catches exceptions because row
    may hide the useful traceback and report only a numeric exit code. On any
    failure, the wrapper writes ``export/export_pose_error.log`` inside the
    action workspace when possible and also prints the same report to stderr.
    This makes direct debugging and grubicy/row debugging use the same error
    information.
    """
    directory: str | Path | None = None

    try:
        if len(argv) != 2:
            raise ValueError(
                "Usage error: expected exactly one argument: the grubicy job directory.\n"
                f"Received argv: {argv!r}"
            )

        directory = argv[1]
        main(directory)
        return 0

    except Exception as exc:  # noqa: BLE001
        report = _build_failure_report(directory=directory, exc=exc)
        log_path = _error_log_path(directory)

        try:
            _write_text(log_path, report)
        except Exception as log_exc:  # noqa: BLE001
            report = (
                report
                + "\n\n"
                + "Could not write export_pose error log\n"
                + "-------------------------------------\n"
                + f"Target log path: {log_path}\n"
                + f"Logging exception: {type(log_exc).__name__}: {log_exc}\n"
            )

        print(report, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli(sys.argv))
