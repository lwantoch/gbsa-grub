"""Run one MD stage from the gbsa-grub grubicy workflow.

This action runner is shared by the MD actions declared in ``pipeline.toml``:
``md_sd``, ``md_cg``, ``md_nvt_res``, ``md_npt_res``, ``md_npt``, and
``md_production``. The current Grubicy action name is read from the job
statepoint and selects the parent manifest, output directory, parameter block,
and gbsa-pipeline MD helper. Each stage loads a parent ``.gro``/``.top`` pair,
runs exactly one MD step, saves a new ``system.gro``/``system.top`` pair, and
writes a small ``result.json`` manifest for the next stage.
"""

from __future__ import annotations

import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from grubicy import open_job_from_directory, parent_file

from config.md_params import (
    CG_PARAMS,
    HEATING_PARAMS,
    NPT_PARAMS,
    NPT_RESTRAINED_PARAMS,
    PRODUCTION_PARAMS,
    SD_PARAMS,
)

BSS: Any = None
run_minimization: Any = None
run_heating: Any = None
run_npt_equilibration: Any = None
run_production: Any = None
load_bss_system_from_gromacs: Any = None
save_bss_system_to_gromacs: Any = None

ERROR_LOG_NAME = "md_stage_error.log"
RUN_LOG_NAME = "md_stage.log"

STAGE_CONFIGS: dict[str, dict[str, Any]] = {
    "md_sd": {
        "parent_manifest": "setup/result.json",
        "output_dir": "sd",
        "runner": "minimization",
        "params": SD_PARAMS,
        "time_ps": None,
        "restraint": None,
        "description": "steepest-descent minimization",
    },
    "md_cg": {
        "parent_manifest": "sd/result.json",
        "output_dir": "cg",
        "runner": "minimization",
        "params": CG_PARAMS,
        "time_ps": None,
        "restraint": None,
        "description": "conjugate-gradient minimization",
    },
    "md_nvt_res": {
        "parent_manifest": "cg/result.json",
        "output_dir": "nvt_res",
        "runner": "heating",
        "params": HEATING_PARAMS,
        "time_ps": 50,
        "restraint": "backbone",
        "description": "restrained NVT heating",
    },
    "md_npt_res": {
        "parent_manifest": "nvt_res/result.json",
        "output_dir": "npt_res",
        "runner": "npt",
        "params": NPT_RESTRAINED_PARAMS,
        "time_ps": 100,
        "restraint": "backbone",
        "description": "restrained NPT equilibration",
    },
    "md_npt": {
        "parent_manifest": "npt_res/result.json",
        "output_dir": "npt",
        "runner": "npt",
        "params": NPT_PARAMS,
        "time_ps": 100,
        "restraint": None,
        "description": "unrestrained NPT equilibration",
    },
    "md_production": {
        "parent_manifest": "npt/result.json",
        "output_dir": "production",
        "runner": "production",
        "params": PRODUCTION_PARAMS,
        "time_ps": 500,
        "restraint": None,
        "description": "production MD",
    },
}


def _print_step(message: str) -> None:
    """Print one plain progress line immediately."""
    print(f"[md_stage] {message}", flush=True)


def _append_text(path: Path, text: str) -> None:
    """Append one line to a local action log."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(text.rstrip() + "\n")


def _log_step(log_path: Path | None, message: str) -> None:
    """Print a progress line and append it to the run log when available."""
    _print_step(message)

    if log_path is not None:
        timestamp = datetime.now(timezone.utc).isoformat()
        _append_text(log_path, f"{timestamp} {message}")


def _read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON object from a stage manifest."""
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
    """Write a deterministic JSON stage manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _load_bss_dependency() -> Any:
    """Import BioSimSpace lazily."""
    global BSS

    if BSS is None:
        import BioSimSpace as imported_bss  # noqa: PLC0415

        BSS = imported_bss

    return BSS


def _load_md_runners() -> None:
    """Import gbsa-pipeline MD runner functions lazily."""
    global run_heating
    global run_minimization
    global run_npt_equilibration
    global run_production

    if (
        run_minimization is not None
        and run_heating is not None
        and run_npt_equilibration is not None
        and run_production is not None
    ):
        return

    from gbsa_pipeline.md import (  # noqa: PLC0415
        run_heating as imported_run_heating,
    )
    from gbsa_pipeline.md import (
        run_minimization as imported_run_minimization,
    )
    from gbsa_pipeline.md import (
        run_npt_equilibration as imported_run_npt_equilibration,
    )
    from gbsa_pipeline.md import (
        run_production as imported_run_production,
    )

    run_minimization = imported_run_minimization
    run_heating = imported_run_heating
    run_npt_equilibration = imported_run_npt_equilibration
    run_production = imported_run_production


def _load_md_io_helpers() -> None:
    """Import gbsa-pipeline GROMACS/BioSimSpace IO helpers lazily."""
    global load_bss_system_from_gromacs
    global save_bss_system_to_gromacs

    if (
        load_bss_system_from_gromacs is not None
        and save_bss_system_to_gromacs is not None
    ):
        return

    from gbsa_pipeline.md_io import (  # noqa: PLC0415
        load_bss_system_from_gromacs as imported_load_bss_system_from_gromacs,
    )
    from gbsa_pipeline.md_io import (
        save_bss_system_to_gromacs as imported_save_bss_system_to_gromacs,
    )

    load_bss_system_from_gromacs = imported_load_bss_system_from_gromacs
    save_bss_system_to_gromacs = imported_save_bss_system_to_gromacs


def _read_parent_manifest(job: Any, relative_path: str) -> dict[str, Any]:
    """Read a stage manifest from the direct Grubicy parent job."""
    manifest_path = parent_file(job, relative_path, must_exist=True)
    return _read_json(manifest_path)


def _required_path(data: dict[str, Any], key: str) -> Path:
    """Return a required manifest path as a resolved file path."""
    if key not in data:
        available = ", ".join(sorted(data)) or "<none>"
        raise KeyError(
            f"Parent manifest is missing required key: {key}\n"
            f"Available keys: {available}\n"
            'Each MD parent manifest must contain "gro_file" and "top_file".'
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


def _current_action(job: Any) -> str:
    """Return the current Grubicy action name from the job statepoint."""
    try:
        action = str(job.sp["action"])
    except KeyError as exc:
        raise KeyError(
            "Current job statepoint is missing required key: action"
        ) from exc

    if action not in STAGE_CONFIGS:
        supported = ", ".join(sorted(STAGE_CONFIGS))
        raise ValueError(
            f"Unsupported MD action '{action}'. Supported actions: {supported}"
        )

    return action


def _run_phase(
    phase_name: str,
    log_path: Path | None,
    func: Any,
    *args: Any,
    **kwargs: Any,
) -> Any:
    """Run one named phase and attach the phase name to failures."""
    _log_step(log_path, f"START {phase_name}")

    try:
        result = func(*args, **kwargs)
    except BaseException as exc:
        _log_step(log_path, f"FAILED {phase_name}: {type(exc).__name__}: {exc}")
        raise RuntimeError(f"MD stage phase failed: {phase_name}") from exc

    _log_step(log_path, f"DONE {phase_name}")
    return result


def _run_stage(
    *,
    action: str,
    system: Any,
    work_dir: Path,
    log_path: Path | None,
) -> Any:
    """Run the gbsa-pipeline MD helper selected by the current action."""
    _load_md_runners()
    bss = _load_bss_dependency()

    stage = STAGE_CONFIGS[action]
    runner = stage["runner"]

    if runner == "minimization":
        return _run_phase(
            f"run {stage['description']}",
            log_path,
            run_minimization,
            system,
            work_dir=work_dir,
            params=stage["params"],
        )

    if runner == "heating":
        return _run_phase(
            f"run {stage['description']}",
            log_path,
            run_heating,
            stage["time_ps"] * bss.Units.Time.picosecond,
            system,
            work_dir=work_dir,
            params=stage["params"],
            temperature_start=50 * bss.Units.Temperature.kelvin,
            temperature_end=300 * bss.Units.Temperature.kelvin,
            restraint=stage["restraint"],
        )

    if runner == "npt":
        return _run_phase(
            f"run {stage['description']}",
            log_path,
            run_npt_equilibration,
            stage["time_ps"] * bss.Units.Time.picosecond,
            system,
            work_dir=work_dir,
            params=stage["params"],
            restraint=stage["restraint"],
        )

    if runner == "production":
        return _run_phase(
            f"run {stage['description']}",
            log_path,
            run_production,
            stage["time_ps"] * bss.Units.Time.picosecond,
            system,
            work_dir=work_dir,
            params=stage["params"],
        )

    raise ValueError(f"Unsupported MD runner for action '{action}': {runner}")


def main(directory: str | Path) -> None:
    """Load the parent system, run one MD stage, and write stage products."""
    action_dir = Path(directory).resolve()
    root_log = action_dir / RUN_LOG_NAME

    _log_step(root_log, f"action directory = {action_dir}")

    job = _run_phase(
        "open grubicy job",
        root_log,
        open_job_from_directory,
        str(action_dir),
    )
    action = _run_phase("read grubicy action name", root_log, _current_action, job)

    stage = STAGE_CONFIGS[action]
    output_dir = Path(job.fn(stage["output_dir"]))
    process_dir = output_dir / "process"
    stage_log = output_dir / RUN_LOG_NAME

    output_dir.mkdir(parents=True, exist_ok=True)
    process_dir.mkdir(parents=True, exist_ok=True)

    _log_step(stage_log, f"action = {action}")
    _log_step(stage_log, f"stage = {stage['description']}")
    _log_step(stage_log, f"parent manifest = {stage['parent_manifest']}")
    _log_step(stage_log, f"output directory = {output_dir}")
    _log_step(stage_log, f"process directory = {process_dir}")

    parent = _run_phase(
        "read parent manifest",
        stage_log,
        _read_parent_manifest,
        job,
        stage["parent_manifest"],
    )
    parent_gro = _run_phase(
        "resolve parent gro_file",
        stage_log,
        _required_path,
        parent,
        "gro_file",
    )
    parent_top = _run_phase(
        "resolve parent top_file",
        stage_log,
        _required_path,
        parent,
        "top_file",
    )

    _log_step(stage_log, f"parent gro_file = {parent_gro}")
    _log_step(stage_log, f"parent top_file = {parent_top}")

    _run_phase("import MD IO helpers", stage_log, _load_md_io_helpers)

    system = _run_phase(
        "load parent GROMACS system",
        stage_log,
        load_bss_system_from_gromacs,
        parent_gro,
        parent_top,
    )

    result_system = _run_stage(
        action=action,
        system=system,
        work_dir=process_dir,
        log_path=stage_log,
    )

    gro_file, top_file = _run_phase(
        "save child GROMACS system",
        stage_log,
        save_bss_system_to_gromacs,
        result_system,
        output_dir / "system",
    )

    result_json = output_dir / "result.json"
    _run_phase(
        "write MD stage manifest",
        stage_log,
        _write_json,
        result_json,
        {
            "action": action,
            "stage": stage["output_dir"],
            "description": stage["description"],
            "parent_gro_file": str(parent_gro),
            "parent_top_file": str(parent_top),
            "gro_file": str(Path(gro_file).resolve()),
            "top_file": str(Path(top_file).resolve()),
            "process_dir": str(process_dir),
            "output_dir": str(output_dir),
            "runner": stage["runner"],
            "time_ps": stage["time_ps"],
            "restraint": stage["restraint"],
        },
    )

    _log_step(stage_log, f"finished action = {action}")
    _log_step(stage_log, f"wrote gro_file = {Path(gro_file).resolve()}")
    _log_step(stage_log, f"wrote top_file = {Path(top_file).resolve()}")
    _log_step(stage_log, f"wrote result_json = {result_json}")


def _error_log_path(directory: str | Path | None) -> Path:
    """Return the best available compact error-log path."""
    if directory is None:
        return Path(ERROR_LOG_NAME).resolve()

    return Path(directory).resolve() / ERROR_LOG_NAME


def _cli(argv: list[str]) -> int:
    """Run the action and write a compact traceback log on failure."""
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

    except BaseException as exc:
        report = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        log_path = _error_log_path(directory)

        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(report, encoding="utf-8")
        except BaseException as log_exc:
            report += (
                "\nCould not write md_stage error log\n"
                f"Target log path: {log_path}\n"
                f"Logging exception: {type(log_exc).__name__}: {log_exc}\n"
            )

        print(report, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli(sys.argv))
