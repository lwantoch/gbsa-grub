"""Unit tests for the setup-system grubicy action.

These tests cover the thin workflow runner that consumes ``export/result.json``
and writes stable setup-stage products for the downstream MD chain. They do not
run OpenMM parametrization, solvation, BioSimSpace, or GROMACS because those
behaviors belong to ``gbsa-pipeline`` library and integration tests. The Grubicy
helper functions and library calls are monkeypatched so the test focuses on
gbsa-grub orchestration wiring. The action intentionally combines
parametrization and solvation in one process because the current solvation
helper needs in-memory objects from parametrization.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from actions import setup_system


class FakeJob:
    """Minimal signac-job stand-in for setup-system action tests.

    The setup action needs ``job.fn(...)`` for output paths and the job object
    for Grubicy helper calls. Runtime action parameters are provided by
    monkeypatching ``load_action_params`` directly, so the fake does not need a
    realistic Grubicy/signac statepoint. This keeps the test independent from an
    initialized project and a materialized workflow. The workspace is backed by
    ``tmp_path`` and is isolated for each test.
    """

    def __init__(self, workspace: Path) -> None:
        self.workspace = workspace

    def fn(self, relative_path: str) -> str:
        """Return a path inside the fake job workspace."""
        return str(self.workspace / relative_path)


def _setup_params() -> setup_system.SetupSystemParams:
    """Return validated setup parameters used by setup-system tests.

    The production action obtains these values through Grubicy's
    ``load_action_params`` helper. Unit tests patch that helper and return this
    model instance directly so they can exercise setup orchestration without a
    real Grubicy workflow. The values match the smoke-workflow solvation
    settings used in the existing tests. Returning the real Pydantic model keeps
    the tests aligned with the action boundary.
    """
    return setup_system.SetupSystemParams(
        water_model="tip3p",
        box_shape="truncated_octahedron",
        padding_nm=1.0,
        ion_concentration=0.15,
        neutralize=True,
    )


def test_setup_system_parametrizes_solvates_and_writes_setup_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the setup action with mocked parametrization and solvation.

    The parent export manifest provides the original protein PDB and the docked
    SDF produced by the pose-export action. The heavy library calls are replaced
    with small fakes that record their arguments and write placeholder setup
    products. The test verifies that ``setup/solvated.gro``,
    ``setup/solvated.top``, and ``setup/result.json`` are created. It also
    checks that the manifest contains the stable ``gro_file`` and ``top_file``
    fields needed by downstream MD actions.
    """
    parent_workspace = tmp_path / "workspace" / "export-job"
    current_workspace = tmp_path / "workspace" / "setup-job"
    export_dir = parent_workspace / "export"
    export_dir.mkdir(parents=True)

    protein_pdb = tmp_path / "inputs" / "protein.pdb"
    docked_sdf = export_dir / "dockligand_vina_out.sdf"
    protein_pdb.parent.mkdir(parents=True)
    protein_pdb.write_text("protein pdb\n", encoding="utf-8")
    docked_sdf.write_text("docked sdf\n", encoding="utf-8")

    parent_manifest_path = export_dir / "result.json"
    parent_manifest_path.write_text(
        json.dumps(
            {
                "protein_pdb": str(protein_pdb),
                "docked_sdf": str(docked_sdf),
                "score": -7.5,
                "rank": 1,
                "engine": "vina",
            }
        ),
        encoding="utf-8",
    )

    job = FakeJob(current_workspace)
    calls: dict[str, Any] = {}

    monkeypatch.setattr(setup_system, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        setup_system,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )
    monkeypatch.setattr(
        setup_system,
        "load_action_params",
        lambda current_job, model: _setup_params(),
    )

    def fake_parametrize(inp: Any) -> Any:
        calls["parametrize"] = inp

        parametrized_gro = (
            current_workspace / "setup" / "parametrization" / "complex.gro"
        )
        parametrized_top = (
            current_workspace / "setup" / "parametrization" / "complex.top"
        )
        crystal_waters = (
            current_workspace / "setup" / "parametrization" / "crystal_waters.pdb"
        )

        parametrized_gro.parent.mkdir(parents=True, exist_ok=True)
        parametrized_gro.write_text("parametrized gro\n", encoding="utf-8")
        parametrized_top.write_text("parametrized top\n", encoding="utf-8")
        crystal_waters.write_text("waters\n", encoding="utf-8")

        return SimpleNamespace(
            gro_file=parametrized_gro,
            top_file=parametrized_top,
            crystal_waters_pdb=crystal_waters,
        )

    def fake_solvate_openmm(
        *,
        parametrized: Any,
        params: Any,
        output_gro: Path,
        output_top: Path,
    ) -> Any:
        calls["solvate_openmm"] = {
            "parametrized": parametrized,
            "params": params,
            "output_gro": output_gro,
            "output_top": output_top,
        }

        output_gro.write_text("solvated gro\n", encoding="utf-8")
        output_top.write_text("solvated top\n", encoding="utf-8")

        return SimpleNamespace(
            gro_file=output_gro,
            top_file=output_top,
        )

    monkeypatch.setattr(setup_system, "parametrize", fake_parametrize)
    monkeypatch.setattr(setup_system, "solvate_openmm", fake_solvate_openmm)

    setup_system.main(current_workspace)

    setup_dir = current_workspace / "setup"
    solvated_gro = setup_dir / "solvated.gro"
    solvated_top = setup_dir / "solvated.top"
    result_json = setup_dir / "result.json"

    assert solvated_gro.exists()
    assert solvated_top.exists()
    assert result_json.exists()

    parametrization_input = calls["parametrize"]

    assert parametrization_input.protein_pdb == protein_pdb.resolve()
    assert parametrization_input.ligand_sdf == docked_sdf.resolve()
    assert parametrization_input.work_dir == setup_dir / "parametrization"

    solvation_call = calls["solvate_openmm"]

    assert solvation_call["output_gro"] == solvated_gro
    assert solvation_call["output_top"] == solvated_top
    assert solvation_call["params"].water_model == "tip3p"
    assert solvation_call["params"].shape == "truncated_octahedron"
    assert solvation_call["params"].padding == 1.0
    assert solvation_call["params"].ion_concentration == 0.15
    assert solvation_call["params"].neutralize is True

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest["protein_pdb"] == str(protein_pdb.resolve())
    assert manifest["docked_sdf"] == str(docked_sdf.resolve())
    assert manifest["solvated_gro"] == str(solvated_gro.resolve())
    assert manifest["solvated_top"] == str(solvated_top.resolve())
    assert manifest["gro_file"] == str(solvated_gro.resolve())
    assert manifest["top_file"] == str(solvated_top.resolve())
    assert manifest["parametrized_gro"] == str(
        (setup_dir / "parametrization" / "complex.gro").resolve()
    )
    assert manifest["parametrized_top"] == str(
        (setup_dir / "parametrization" / "complex.top").resolve()
    )
    assert manifest["crystal_waters_pdb"] == str(
        (setup_dir / "parametrization" / "crystal_waters.pdb").resolve()
    )
    assert manifest["setup_dir"] == str(setup_dir)
    assert manifest["parametrization_dir"] == str(setup_dir / "parametrization")
    assert manifest["solvation_dir"] == str(setup_dir / "solvation")
    assert manifest["water_model"] == "tip3p"
    assert manifest["box_shape"] == "truncated_octahedron"
    assert manifest["padding_nm"] == 1.0
    assert manifest["ion_concentration"] == 0.15
    assert manifest["neutralize"] is True
    assert manifest["score"] == -7.5
    assert manifest["rank"] == 1
    assert manifest["engine"] == "vina"


def test_setup_system_reports_missing_parent_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail early when the export-stage manifest is missing.

    The setup action should not start parametrization when the declared parent
    product is unavailable. This test monkeypatches Grubicy's ``parent_file``
    helper to report that ``export/result.json`` does not exist. The action
    should raise a clear file error before calling any heavy gbsa-pipeline
    function. This keeps missing dependency products visible at the workflow
    boundary.
    """
    job = FakeJob(tmp_path / "workspace" / "setup-job")

    monkeypatch.setattr(setup_system, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        setup_system,
        "load_action_params",
        lambda current_job, model: _setup_params(),
    )

    def missing_parent_file(
        current_job: Any,
        relative_path: str,
        *,
        must_exist: bool = True,
    ) -> Path:
        raise FileNotFoundError(f"Missing parent file: {relative_path}")

    monkeypatch.setattr(setup_system, "parent_file", missing_parent_file)

    with pytest.raises(FileNotFoundError, match=setup_system.PARENT_MANIFEST):
        setup_system.main(job.workspace)


def test_setup_system_reports_missing_required_parent_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail when a file referenced by the parent manifest is missing.

    The setup action needs both the original protein PDB and the exported docked
    SDF before parametrization can begin. This test provides an export manifest
    whose ``docked_sdf`` points to a missing file while the protein file exists.
    The action should fail before parametrization or solvation starts. This
    gives a direct path-level error instead of a later chemistry-tool failure.
    """
    parent_workspace = tmp_path / "workspace" / "export-job"
    current_workspace = tmp_path / "workspace" / "setup-job"
    export_dir = parent_workspace / "export"
    export_dir.mkdir(parents=True)

    protein_pdb = tmp_path / "inputs" / "protein.pdb"
    protein_pdb.parent.mkdir(parents=True)
    protein_pdb.write_text("protein pdb\n", encoding="utf-8")

    parent_manifest_path = export_dir / "result.json"
    parent_manifest_path.write_text(
        json.dumps(
            {
                "protein_pdb": str(protein_pdb),
                "docked_sdf": str(export_dir / "missing.sdf"),
            }
        ),
        encoding="utf-8",
    )

    job = FakeJob(current_workspace)

    monkeypatch.setattr(setup_system, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        setup_system,
        "load_action_params",
        lambda current_job, model: _setup_params(),
    )
    monkeypatch.setattr(
        setup_system,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )

    with pytest.raises(FileNotFoundError, match="docked_sdf"):
        setup_system.main(current_workspace)
