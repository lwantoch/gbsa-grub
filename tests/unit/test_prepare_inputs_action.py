"""Unit tests for the prepare-inputs grubicy action.

These tests cover the thin workflow runner that prepares docking inputs from a
Grubicy action parameter set. They do not run Meeko, receptor preparation,
docking, or any molecular parsing because those behaviors belong to
``gbsa-pipeline`` library tests. The runner is tested by monkeypatching the
library calls and checking that the expected workflow products and manifest
fields are written. This keeps the test focused on gbsa-grub orchestration
wiring.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from actions import prepare_inputs


class FakeJob:
    """Minimal signac-job stand-in for the prepare-inputs action tests.

    The action only needs ``job.fn(...)`` for paths inside the job workspace.
    Runtime action parameters are provided by monkeypatching Grubicy's
    ``load_action_params`` helper directly, so the fake does not need a realistic
    signac statepoint. It intentionally implements no other signac behavior so
    accidental extra runner responsibilities fail during the test. The workspace
    is provided by ``tmp_path`` and therefore isolated for each test case.
    """

    def __init__(self, workspace: Path) -> None:
        self.workspace = workspace

    def fn(self, relative_path: str) -> str:
        """Return a path inside the fake job workspace."""
        return str(self.workspace / relative_path)


def _prepare_params(
    *,
    protein_pdb: Path,
    ligand_sdf: Path,
    box_center: tuple[float, float, float] = (10.115, 39.148, 53.112),
    box_size: tuple[float, float, float] = (10.0, 10.0, 10.0),
) -> prepare_inputs.PrepareInputsParams:
    """Return validated prepare-action parameters for unit tests.

    The production action obtains this model through Grubicy's
    ``load_action_params`` helper. Unit tests patch that helper and return a real
    ``PrepareInputsParams`` instance so they exercise runner orchestration
    without materializing a Grubicy workflow. Input files are still created on
    disk by each test because the action performs explicit file-existence checks.
    The default box values match the smoke-workflow geometry used by the existing
    tests.
    """
    return prepare_inputs.PrepareInputsParams(
        protein_pdb=protein_pdb,
        ligand_sdf=ligand_sdf,
        box_center=box_center,
        box_size=box_size,
    )


def test_prepare_inputs_writes_declared_outputs_and_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the prepare-inputs action with mocked library calls.

    The typed action parameters provide the input files and docking box values
    exactly like a Grubicy root action would. The heavy library functions are
    replaced with small fakes that record their arguments, write placeholder
    products, and return the expected output paths. The test verifies that the
    action creates ``prepare/ligand.pdbqt``, ``prepare/receptor.pdbqt``, and
    ``prepare/result.json``. It also checks that the manifest contains resolved
    input paths and the converted docking box values.
    """
    protein_pdb = tmp_path / "input" / "protein.pdb"
    ligand_sdf = tmp_path / "input" / "ligand.sdf"
    protein_pdb.parent.mkdir(parents=True)
    protein_pdb.write_text("protein\n", encoding="utf-8")
    ligand_sdf.write_text("ligand\n", encoding="utf-8")

    job = FakeJob(workspace=tmp_path / "workspace" / "job-id")
    params = _prepare_params(protein_pdb=protein_pdb, ligand_sdf=ligand_sdf)

    ligand_molecule = object()
    calls: dict[str, Any] = {}

    monkeypatch.setattr(
        prepare_inputs, "open_job_from_directory", lambda directory: job
    )
    monkeypatch.setattr(
        prepare_inputs,
        "load_action_params",
        lambda current_job, model: params,
    )
    monkeypatch.setattr(
        prepare_inputs.shutil,
        "which",
        lambda binary: f"/usr/bin/{binary}",
    )

    def fake_load_first_sdf_molecule(path: Path, *, remove_hs: bool) -> object:
        calls["load_first_sdf_molecule"] = {
            "path": path,
            "remove_hs": remove_hs,
        }
        return ligand_molecule

    def fake_prepare_ligand_with_meeko(
        molecule: object,
        output_path: Path,
        *,
        name: str | None = None,
    ) -> Path:
        calls["prepare_ligand_with_meeko"] = {
            "molecule": molecule,
            "output_path": output_path,
            "name": name,
        }
        output_path.write_text("ligand pdbqt\n", encoding="utf-8")
        return output_path

    def fake_convert_receptor_pdb_to_pdbqt(
        receptor_pdb: Path,
        *,
        output_path: Path,
        mk_prepare_receptor_binary: str,
    ) -> Path:
        calls["convert_receptor_pdb_to_pdbqt"] = {
            "receptor_pdb": receptor_pdb,
            "output_path": output_path,
            "mk_prepare_receptor_binary": mk_prepare_receptor_binary,
        }
        output_path.write_text("receptor pdbqt\n", encoding="utf-8")
        return output_path

    monkeypatch.setattr(
        prepare_inputs,
        "load_first_sdf_molecule",
        fake_load_first_sdf_molecule,
    )
    monkeypatch.setattr(
        prepare_inputs,
        "prepare_ligand_with_meeko",
        fake_prepare_ligand_with_meeko,
    )
    monkeypatch.setattr(
        prepare_inputs,
        "convert_receptor_pdb_to_pdbqt",
        fake_convert_receptor_pdb_to_pdbqt,
    )

    prepare_inputs.main(job.workspace)

    prepare_dir = job.workspace / "prepare"
    ligand_pdbqt = prepare_dir / "ligand.pdbqt"
    receptor_pdbqt = prepare_dir / "receptor.pdbqt"
    result_json = prepare_dir / "result.json"

    assert ligand_pdbqt.exists()
    assert receptor_pdbqt.exists()
    assert result_json.exists()

    assert calls["load_first_sdf_molecule"] == {
        "path": ligand_sdf.resolve(),
        "remove_hs": False,
    }
    assert calls["prepare_ligand_with_meeko"] == {
        "molecule": ligand_molecule,
        "output_path": ligand_pdbqt,
        "name": "DOCKLIG",
    }
    assert calls["convert_receptor_pdb_to_pdbqt"] == {
        "receptor_pdb": protein_pdb.resolve(),
        "output_path": receptor_pdbqt,
        "mk_prepare_receptor_binary": prepare_inputs.MEEKO_RECEPTOR_BINARY,
    }

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest == {
        "box_center": [10.115, 39.148, 53.112],
        "box_size": [10.0, 10.0, 10.0],
        "ligand_pdbqt": str(ligand_pdbqt.resolve()),
        "ligand_sdf": str(ligand_sdf.resolve()),
        "mk_prepare_receptor_binary": prepare_inputs.MEEKO_RECEPTOR_BINARY,
        "prepare_dir": str(prepare_dir),
        "protein_pdb": str(protein_pdb.resolve()),
        "receptor_pdbqt": str(receptor_pdbqt.resolve()),
    }


def test_prepare_inputs_reports_missing_receptor_binary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail before file preparation when the receptor-preparation binary is missing.

    The action uses the receptor preparation binary before writing workflow
    products, so a missing executable should produce a direct runtime error. The
    fake job has valid typed parameters, but the binary check is monkeypatched
    to fail. This keeps the error path independent from the library helpers. The
    message includes the binary name so the workflow failure is actionable.
    """
    protein_pdb = tmp_path / "protein.pdb"
    ligand_sdf = tmp_path / "ligand.sdf"
    protein_pdb.write_text("protein\n", encoding="utf-8")
    ligand_sdf.write_text("ligand\n", encoding="utf-8")

    job = FakeJob(workspace=tmp_path / "workspace" / "job-id")
    params = _prepare_params(
        protein_pdb=protein_pdb,
        ligand_sdf=ligand_sdf,
        box_center=(1.0, 2.0, 3.0),
        box_size=(10.0, 10.0, 10.0),
    )

    monkeypatch.setattr(
        prepare_inputs, "open_job_from_directory", lambda directory: job
    )
    monkeypatch.setattr(
        prepare_inputs,
        "load_action_params",
        lambda current_job, model: params,
    )
    monkeypatch.setattr(prepare_inputs.shutil, "which", lambda binary: None)

    with pytest.raises(RuntimeError, match=prepare_inputs.MEEKO_RECEPTOR_BINARY):
        prepare_inputs.main(job.workspace)


def test_prepare_inputs_rejects_invalid_box_values(tmp_path: Path) -> None:
    """Reject malformed docking box action-parameter values.

    The docking box center and size must each contain exactly three numeric
    values before they are passed to ``DockingBox``. The old local conversion
    helper was removed when the action switched to Grubicy typed parameters, so
    this test now checks the Pydantic parameter model directly. Failing here
    gives a clearer error than letting a later docking step receive malformed
    geometry. The test does not need to call ``main`` because malformed action
    parameters are rejected before runner orchestration starts.
    """
    protein_pdb = tmp_path / "protein.pdb"
    ligand_sdf = tmp_path / "ligand.sdf"

    with pytest.raises(ValidationError):
        prepare_inputs.PrepareInputsParams(
            protein_pdb=protein_pdb,
            ligand_sdf=ligand_sdf,
            box_center=(1.0, 2.0),
            box_size=(10.0, 10.0, 10.0),
        )
