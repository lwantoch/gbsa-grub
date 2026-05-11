# tests/unit/test_md_stage_action.py
"""Unit tests for the shared MD-stage grubicy action.

These tests cover the thin workflow runner that is reused by all MD actions in
the gbsa-grub workflow. They do not run BioSimSpace, GROMACS, minimization,
equilibration, or production MD because those behaviors belong to gbsa-pipeline
library and integration tests. The Grubicy helper functions, system IO helpers,
and MD runner functions are monkeypatched so the tests focus on orchestration
wiring. The important contract is that each action loads the parent gro/top
system, runs exactly one configured MD stage, writes a new system gro/top pair,
and produces a result manifest for the next action.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from actions import md_stage


class FakeJob:
    """Minimal signac-job stand-in for MD-stage action tests.

    The shared MD runner needs ``job.sp`` for the current Grubicy action name,
    ``job.fn(...)`` for output paths, and the job object for Grubicy parent-file
    helper calls. This fake avoids materializing a real signac or Grubicy
    project while preserving the small interface used by the runner. It
    intentionally implements no other behavior so accidental extra runner
    responsibilities become visible during tests. The workspace is backed by
    ``tmp_path`` and is isolated for each test.
    """

    def __init__(self, workspace: Path, action: str) -> None:
        self.workspace = workspace
        self.sp = {"action": action}

    def fn(self, relative_path: str) -> str:
        """Return a path inside the fake job workspace."""
        return str(self.workspace / relative_path)


class FakeBSS:
    """Tiny BioSimSpace unit stand-in for MD-stage dispatch tests.

    The shared runner multiplies stage durations by BioSimSpace unit constants.
    The tests do not need real unit objects, only values that support numeric
    multiplication and can be passed through to fake runner functions. Using
    ``1`` keeps the resulting values equal to the configured picoseconds and
    kelvin values. This avoids importing BioSimSpace while still exercising the
    action dispatch code path.
    """

    class Units:
        """Minimal unit namespace used by the runner."""

        class Time:
            """Minimal time unit namespace."""

            picosecond = 1

        class Temperature:
            """Minimal temperature unit namespace."""

            kelvin = 1


def _write_parent_manifest(
    parent_workspace: Path,
    manifest_path: str,
    *,
    gro_file: Path,
    top_file: Path,
) -> Path:
    """Write a parent stage manifest for an MD action test.

    The shared MD action reads the direct parent manifest declared in its stage
    config. This helper creates that file under the fake parent workspace and
    stores the parent ``gro_file`` and ``top_file`` paths expected by the action.
    It keeps each test focused on runner behavior instead of repetitive JSON
    setup. The referenced files must be created separately so the action's path
    validation remains meaningful. The created manifest path is returned so the
    fake Grubicy ``parent_file`` helper can resolve it directly.
    """
    path = parent_workspace / manifest_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "gro_file": str(gro_file),
                "top_file": str(top_file),
            }
        ),
        encoding="utf-8",
    )
    return path


def test_md_stage_runs_sd_stage_and_writes_declared_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the ``md_sd`` action with mocked IO and MD execution.

    The parent setup manifest provides the solvated ``.gro`` and ``.top`` files.
    The BioSimSpace system loader, minimization runner, and system saver are
    replaced with small fakes that record their arguments and write placeholder
    products. The test verifies that the runner creates ``sd/system.gro``,
    ``sd/system.top``, and ``sd/result.json``. It also checks that the manifest
    records the parent files, output files, process directory, action name, and
    selected stage runner.
    """
    parent_workspace = tmp_path / "workspace" / "setup-job"
    current_workspace = tmp_path / "workspace" / "md-sd-job"

    parent_gro = parent_workspace / "setup" / "solvated.gro"
    parent_top = parent_workspace / "setup" / "solvated.top"
    parent_gro.parent.mkdir(parents=True)
    parent_gro.write_text("parent gro\n", encoding="utf-8")
    parent_top.write_text("parent top\n", encoding="utf-8")

    parent_manifest_path = _write_parent_manifest(
        parent_workspace,
        "setup/result.json",
        gro_file=parent_gro,
        top_file=parent_top,
    )

    job = FakeJob(current_workspace, action="md_sd")
    loaded_system = object()
    minimized_system = object()
    calls: dict[str, Any] = {}

    monkeypatch.setattr(md_stage, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        md_stage,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )
    monkeypatch.setattr(md_stage, "_load_md_io_helpers", lambda: None)
    monkeypatch.setattr(md_stage, "_load_md_runners", lambda: None)
    monkeypatch.setattr(md_stage, "_load_bss_dependency", lambda: FakeBSS)

    def fake_load_bss_system_from_gromacs(gro_file: Path, top_file: Path) -> object:
        calls["load"] = {
            "gro_file": gro_file,
            "top_file": top_file,
        }
        return loaded_system

    def fake_run_minimization(
        system: object,
        *,
        work_dir: Path,
        params: dict[str, Any],
    ) -> object:
        calls["run_minimization"] = {
            "system": system,
            "work_dir": work_dir,
            "params": params,
        }
        return minimized_system

    def fake_save_bss_system_to_gromacs(
        system: object,
        output_prefix: Path,
    ) -> tuple[Path, Path]:
        calls["save"] = {
            "system": system,
            "output_prefix": output_prefix,
        }
        gro_file = output_prefix.with_suffix(".gro")
        top_file = output_prefix.with_suffix(".top")
        gro_file.parent.mkdir(parents=True, exist_ok=True)
        gro_file.write_text("saved gro\n", encoding="utf-8")
        top_file.write_text("saved top\n", encoding="utf-8")
        return gro_file, top_file

    monkeypatch.setattr(
        md_stage,
        "load_bss_system_from_gromacs",
        fake_load_bss_system_from_gromacs,
    )
    monkeypatch.setattr(md_stage, "run_minimization", fake_run_minimization)
    monkeypatch.setattr(
        md_stage,
        "save_bss_system_to_gromacs",
        fake_save_bss_system_to_gromacs,
    )

    md_stage.main(current_workspace)

    output_dir = current_workspace / "sd"
    process_dir = output_dir / "process"
    system_gro = output_dir / "system.gro"
    system_top = output_dir / "system.top"
    result_json = output_dir / "result.json"

    assert system_gro.exists()
    assert system_top.exists()
    assert result_json.exists()

    assert calls["load"] == {
        "gro_file": parent_gro.resolve(),
        "top_file": parent_top.resolve(),
    }
    assert calls["run_minimization"] == {
        "system": loaded_system,
        "work_dir": process_dir,
        "params": md_stage.SD_PARAMS,
    }
    assert calls["save"] == {
        "system": minimized_system,
        "output_prefix": output_dir / "system",
    }

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest == {
        "action": "md_sd",
        "description": "steepest-descent minimization",
        "gro_file": str(system_gro.resolve()),
        "output_dir": str(output_dir),
        "parent_gro_file": str(parent_gro.resolve()),
        "parent_top_file": str(parent_top.resolve()),
        "process_dir": str(process_dir),
        "restraint": None,
        "runner": "minimization",
        "stage": "sd",
        "time_ps": None,
        "top_file": str(system_top.resolve()),
    }


@pytest.mark.parametrize(
    ("action", "parent_manifest", "output_dir_name", "runner_name"),
    [
        ("md_cg", "sd/result.json", "cg", "minimization"),
        ("md_nvt_res", "cg/result.json", "nvt_res", "heating"),
        ("md_npt_res", "nvt_res/result.json", "npt_res", "npt"),
        ("md_npt", "npt_res/result.json", "npt", "npt"),
        ("md_production", "npt/result.json", "production", "production"),
    ],
)
def test_md_stage_dispatches_supported_actions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    action: str,
    parent_manifest: str,
    output_dir_name: str,
    runner_name: str,
) -> None:
    """Check that supported MD actions dispatch to the expected runner.

    The shared runner uses ``job.sp["action"]`` to select the stage-specific
    parent manifest, output directory, parameter block, and gbsa-pipeline MD
    helper. This test exercises the non-SD actions with mocked MD functions so
    no simulation is started. It verifies that each action writes the expected
    output directory and reports the selected runner in its result manifest.
    Detailed parameter correctness is covered by the stage configuration itself
    and by library-level MD tests.
    """
    parent_workspace = tmp_path / "workspace" / "parent-job"
    current_workspace = tmp_path / "workspace" / f"{action}-job"

    parent_gro = parent_workspace / "parent" / "system.gro"
    parent_top = parent_workspace / "parent" / "system.top"
    parent_gro.parent.mkdir(parents=True)
    parent_gro.write_text("parent gro\n", encoding="utf-8")
    parent_top.write_text("parent top\n", encoding="utf-8")

    parent_manifest_path = _write_parent_manifest(
        parent_workspace,
        parent_manifest,
        gro_file=parent_gro,
        top_file=parent_top,
    )

    job = FakeJob(current_workspace, action=action)
    loaded_system = object()
    result_system = object()
    calls: dict[str, Any] = {}

    monkeypatch.setattr(md_stage, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        md_stage,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )
    monkeypatch.setattr(md_stage, "_load_md_io_helpers", lambda: None)
    monkeypatch.setattr(md_stage, "_load_md_runners", lambda: None)
    monkeypatch.setattr(md_stage, "_load_bss_dependency", lambda: FakeBSS)
    monkeypatch.setattr(
        md_stage,
        "load_bss_system_from_gromacs",
        lambda gro_file, top_file: loaded_system,
    )

    def fake_stage_runner(*args: Any, **kwargs: Any) -> object:
        calls["runner"] = {
            "args": args,
            "kwargs": kwargs,
        }
        return result_system

    if runner_name == "minimization":
        monkeypatch.setattr(md_stage, "run_minimization", fake_stage_runner)
    elif runner_name == "heating":
        monkeypatch.setattr(md_stage, "run_heating", fake_stage_runner)
    elif runner_name == "npt":
        monkeypatch.setattr(md_stage, "run_npt_equilibration", fake_stage_runner)
    elif runner_name == "production":
        monkeypatch.setattr(md_stage, "run_production", fake_stage_runner)
    else:  # pragma: no cover - defensive against typo in parametrization
        raise AssertionError(f"Unsupported runner name in test: {runner_name}")

    def fake_save_bss_system_to_gromacs(
        system: object,
        output_prefix: Path,
    ) -> tuple[Path, Path]:
        calls["saved_system"] = system
        gro_file = output_prefix.with_suffix(".gro")
        top_file = output_prefix.with_suffix(".top")
        gro_file.parent.mkdir(parents=True, exist_ok=True)
        gro_file.write_text("saved gro\n", encoding="utf-8")
        top_file.write_text("saved top\n", encoding="utf-8")
        return gro_file, top_file

    monkeypatch.setattr(
        md_stage,
        "save_bss_system_to_gromacs",
        fake_save_bss_system_to_gromacs,
    )

    md_stage.main(current_workspace)

    output_dir = current_workspace / output_dir_name
    result_json = output_dir / "result.json"

    assert (output_dir / "system.gro").exists()
    assert (output_dir / "system.top").exists()
    assert result_json.exists()
    assert calls["saved_system"] is result_system

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest["action"] == action
    assert manifest["stage"] == output_dir_name
    assert manifest["runner"] == runner_name
    assert manifest["gro_file"] == str((output_dir / "system.gro").resolve())
    assert manifest["top_file"] == str((output_dir / "system.top").resolve())
    assert "runner" in calls


def test_md_stage_reports_missing_parent_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail early when the expected direct parent manifest is missing.

    A child MD action should not try to load a system when the declared parent
    product is unavailable. This test monkeypatches Grubicy's ``parent_file``
    helper to report that the parent manifest does not exist. The action should
    fail before calling the BioSimSpace loader or any MD runner. The runner wraps
    phase failures, so the original file error remains available as the chained
    cause.
    """
    job = FakeJob(tmp_path / "workspace" / "md-job", action="md_sd")

    monkeypatch.setattr(md_stage, "open_job_from_directory", lambda directory: job)

    def missing_parent_file(
        current_job: Any,
        relative_path: str,
        *,
        must_exist: bool = True,
    ) -> Path:
        raise FileNotFoundError(f"Missing parent file: {relative_path}")

    monkeypatch.setattr(md_stage, "parent_file", missing_parent_file)

    with pytest.raises(RuntimeError, match="read parent manifest") as exc_info:
        md_stage.main(job.workspace)

    assert isinstance(exc_info.value.__cause__, FileNotFoundError)
    assert "setup/result.json" in str(exc_info.value.__cause__)


def test_md_stage_reports_unsupported_action(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail before IO when the Grubicy action name is unsupported.

    The shared runner relies on ``job.sp["action"]`` to choose one of the known
    MD stages. If a job was materialized with an unexpected action name, the
    runner should fail before reading parent files or launching MD. This makes
    pipeline specification mistakes easier to diagnose. The runner wraps phase
    failures, so the original value error remains available as the chained
    cause.
    """
    job = FakeJob(tmp_path / "workspace" / "bad-job", action="md_unknown")

    monkeypatch.setattr(md_stage, "open_job_from_directory", lambda directory: job)

    with pytest.raises(RuntimeError, match="read grubicy action name") as exc_info:
        md_stage.main(job.workspace)

    assert isinstance(exc_info.value.__cause__, ValueError)
    assert "Unsupported MD action" in str(exc_info.value.__cause__)
