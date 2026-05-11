"""Unit tests for the docking grubicy action.

These tests cover the thin workflow runner that consumes ``prepare/result.json``
and writes docking-stage products. They do not run AutoDock Vina and do not
materialize a real grubicy workflow, because the behavior under test is runner
wiring rather than docking science. The Vina engine and grubicy helper functions
are monkeypatched so the test can verify parent-manifest reading, output
normalization, and result-manifest contents. Library-level docking behavior
belongs to gbsa-pipeline tests.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from actions import dock


class FakeJob:
    """Minimal signac-job stand-in for the docking action tests.

    The docking action only needs ``job.fn(...)`` and access to the current job
    object for Grubicy parent-file helper calls. This fake keeps the test
    independent from an initialized signac project and a materialized Grubicy
    workflow. It intentionally implements only the small interface used by the
    action so new runner responsibilities become visible in the test. The
    workspace is backed by ``tmp_path`` and is therefore isolated for each test.
    """

    def __init__(self, workspace: Path) -> None:
        self.workspace = workspace

    def fn(self, relative_path: str) -> str:
        """Return a path inside the fake job workspace."""
        return str(self.workspace / relative_path)


def test_dock_action_reads_parent_manifest_runs_engine_and_writes_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the docking action with mocked Grubicy and Vina dependencies.

    The parent prepare manifest is written exactly where the fake Grubicy
    ``parent_file`` helper returns it. The fake Vina engine writes the library's
    ligand-stem-based output names, and the action then normalizes them to the
    declared ``dockligand_*`` workflow product names. The test verifies that the
    final pose, final log, and ``docking/result.json`` exist and contain the
    expected paths and score metadata. No real docking binary is executed.
    """
    parent_workspace = tmp_path / "workspace" / "parent-job"
    current_workspace = tmp_path / "workspace" / "dock-job"
    prepare_dir = parent_workspace / "prepare"
    prepare_dir.mkdir(parents=True)

    protein_pdb = prepare_dir / "protein.pdb"
    ligand_sdf = prepare_dir / "ligand.sdf"
    ligand_pdbqt = prepare_dir / "ligand.pdbqt"
    receptor_pdbqt = prepare_dir / "receptor.pdbqt"

    protein_pdb.write_text("protein pdb\n", encoding="utf-8")
    ligand_sdf.write_text("ligand sdf\n", encoding="utf-8")
    ligand_pdbqt.write_text("ligand pdbqt\n", encoding="utf-8")
    receptor_pdbqt.write_text("receptor pdbqt\n", encoding="utf-8")

    parent_manifest = {
        "protein_pdb": str(protein_pdb),
        "ligand_sdf": str(ligand_sdf),
        "ligand_pdbqt": str(ligand_pdbqt),
        "receptor_pdbqt": str(receptor_pdbqt),
        "box_center": [10.115, 39.148, 53.112],
        "box_size": [10.0, 10.0, 10.0],
    }
    parent_manifest_path = prepare_dir / "result.json"
    parent_manifest_path.write_text(json.dumps(parent_manifest), encoding="utf-8")

    job = FakeJob(current_workspace)
    calls: dict[str, Any] = {}

    monkeypatch.setattr(dock, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(dock.shutil, "which", lambda binary: f"/usr/bin/{binary}")
    monkeypatch.setattr(
        dock,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )

    class FakeVinaEngine:
        """Small fake for the library Vina engine used by the action.

        The real Vina engine is covered in gbsa-pipeline tests and should not be
        launched from this runner unit test. This fake records the request object
        so the test can inspect which receptor, ligand, box, and workdir were
        passed by the action. It writes placeholder generated products using the
        same ligand-stem naming pattern as the library engine. The returned
        object exposes only the attributes consumed by the runner.
        """

        def __init__(self, *, binary: str) -> None:
            calls["binary"] = binary

        def dock(self, *, request: Any) -> Any:
            calls["request"] = request

            generated_pose = Path(request.workdir) / "ligand_vina_out.pdbqt"
            generated_log = Path(request.workdir) / "ligand_vina.log"
            generated_pose.write_text("pose pdbqt\n", encoding="utf-8")
            generated_log.write_text("vina log\n", encoding="utf-8")

            pose = SimpleNamespace(
                pose_path=generated_pose,
                score=-7.5,
                rank=1,
                metadata={
                    "returncode": 0,
                    "log_file": str(generated_log),
                },
            )

            return SimpleNamespace(
                poses=[pose],
                engine="vina",
            )

    monkeypatch.setattr(dock, "VinaEngine", FakeVinaEngine)

    dock.main(current_workspace)

    docking_dir = current_workspace / "docking"
    final_pose = docking_dir / "dockligand_vina_out.pdbqt"
    final_log = docking_dir / "dockligand_vina.log"
    result_json = docking_dir / "result.json"

    assert final_pose.exists()
    assert final_log.exists()
    assert result_json.exists()
    assert not (docking_dir / "ligand_vina_out.pdbqt").exists()
    assert not (docking_dir / "ligand_vina.log").exists()

    request = calls["request"]

    assert calls["binary"] == "vina"
    assert request.receptor == receptor_pdbqt.resolve()
    assert request.ligands == [ligand_pdbqt.resolve()]
    assert request.workdir == docking_dir
    assert list(request.box.center) == [10.115, 39.148, 53.112]
    assert list(request.box.size) == [10.0, 10.0, 10.0]

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest == {
        "protein_pdb": str(protein_pdb.resolve()),
        "ligand_sdf": str(ligand_sdf.resolve()),
        "box_center": [10.115, 39.148, 53.112],
        "box_size": [10.0, 10.0, 10.0],
        "docking_dir": str(docking_dir),
        "engine": "vina",
        "ligand_pdbqt": str(ligand_pdbqt.resolve()),
        "pose_pdbqt": str(final_pose.resolve()),
        "rank": 1,
        "receptor_pdbqt": str(receptor_pdbqt.resolve()),
        "score": -7.5,
        "vina_log": str(final_log.resolve()),
    }


def test_dock_action_reports_missing_parent_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail early when the prepare-stage manifest is not available.

    A child Grubicy action should not guess upstream paths when the declared
    parent product is missing. This test monkeypatches Grubicy's ``parent_file``
    helper to report that ``prepare/result.json`` does not exist. The action
    should fail before Vina is constructed or launched. This makes missing
    dependency products easier to diagnose from the workflow layer.
    """
    job = FakeJob(tmp_path / "workspace" / "dock-job")

    monkeypatch.setattr(dock, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(dock.shutil, "which", lambda binary: f"/usr/bin/{binary}")

    def missing_parent_file(
        current_job: Any,
        relative_path: str,
        *,
        must_exist: bool = True,
    ) -> Path:
        raise FileNotFoundError(f"Missing parent file: {relative_path}")

    monkeypatch.setattr(dock, "parent_file", missing_parent_file)

    with pytest.raises(FileNotFoundError, match=dock.PARENT_MANIFEST):
        dock.main(job.workspace)


def test_dock_action_reports_vina_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail when the library Vina engine reports a non-zero return code.

    The runner should not write a success manifest when Vina failed. This test
    provides valid parent files and a fake engine result with ``returncode`` set
    to a non-zero value. The action should raise a runtime error that includes
    enough context to inspect the Vina log. The actual interpretation of Vina
    stderr remains a gbsa-pipeline responsibility.
    """
    parent_workspace = tmp_path / "workspace" / "parent-job"
    current_workspace = tmp_path / "workspace" / "dock-job"
    prepare_dir = parent_workspace / "prepare"
    prepare_dir.mkdir(parents=True)

    protein_pdb = prepare_dir / "protein.pdb"
    ligand_sdf = prepare_dir / "ligand.sdf"
    ligand_pdbqt = prepare_dir / "ligand.pdbqt"
    receptor_pdbqt = prepare_dir / "receptor.pdbqt"

    protein_pdb.write_text("protein pdb\n", encoding="utf-8")
    ligand_sdf.write_text("ligand sdf\n", encoding="utf-8")
    ligand_pdbqt.write_text("ligand pdbqt\n", encoding="utf-8")
    receptor_pdbqt.write_text("receptor pdbqt\n", encoding="utf-8")

    parent_manifest_path = prepare_dir / "result.json"
    parent_manifest_path.write_text(
        json.dumps(
            {
                "protein_pdb": str(protein_pdb),
                "ligand_sdf": str(ligand_sdf),
                "ligand_pdbqt": str(ligand_pdbqt),
                "receptor_pdbqt": str(receptor_pdbqt),
                "box_center": [1.0, 2.0, 3.0],
                "box_size": [10.0, 10.0, 10.0],
            }
        ),
        encoding="utf-8",
    )

    job = FakeJob(current_workspace)

    monkeypatch.setattr(dock, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(dock.shutil, "which", lambda binary: f"/usr/bin/{binary}")
    monkeypatch.setattr(
        dock,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )

    class FailingVinaEngine:
        """Fake Vina engine that reports a failed docking run."""

        def __init__(self, *, binary: str) -> None:
            self.binary = binary

        def dock(self, *, request: Any) -> Any:
            pose = SimpleNamespace(
                pose_path=Path(request.workdir) / "ligand_vina_out.pdbqt",
                score=None,
                rank=None,
                metadata={
                    "returncode": 1,
                    "log_file": str(Path(request.workdir) / "ligand_vina.log"),
                },
            )
            return SimpleNamespace(poses=[pose], engine="vina")

    monkeypatch.setattr(dock, "VinaEngine", FailingVinaEngine)

    with pytest.raises(RuntimeError, match="Vina docking failed"):
        dock.main(current_workspace)

    assert not (current_workspace / "docking" / "result.json").exists()
