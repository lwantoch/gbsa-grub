# tests/unit/test_export_pose_action.py
"""Unit tests for the pose-export grubicy action.

These tests cover the thin workflow runner that consumes ``docking/result.json``
and writes export-stage products. They do not run Meeko, parse real PDBQT files,
or reconstruct ligand chemistry because those behaviors belong to
``gbsa-pipeline`` tests. The Grubicy helper functions and library conversion
calls are monkeypatched so the test can focus on manifest wiring, stable output
names, and failure handling. This keeps gbsa-grub tests focused on orchestration
rather than scientific file conversion.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from actions import export_pose


class FakeJob:
    """Minimal signac-job stand-in for export action tests.

    The export action only needs ``job.fn(...)`` and the job object for Grubicy
    parent-file helper calls. This fake avoids materializing a real signac or
    Grubicy project while preserving the small interface the runner uses. It
    intentionally implements no other behavior so accidental new runner
    responsibilities become visible in tests. The workspace is backed by
    ``tmp_path`` and is isolated for each test.
    """

    def __init__(self, workspace: Path) -> None:
        self.workspace = workspace

    def fn(self, relative_path: str) -> str:
        """Return a path inside the fake job workspace."""
        return str(self.workspace / relative_path)


def test_export_pose_reads_parent_manifest_exports_sdf_and_writes_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the export action with mocked Grubicy and conversion helpers.

    The parent docking manifest contains the docked PDBQT path, the original
    ligand SDF template path, and score metadata from Vina. The heavy library
    calls are replaced with small fakes that record their arguments and write a
    placeholder SDF product. The test verifies that the declared
    ``export/dockligand_vina_out.sdf`` and ``export/result.json`` products are
    created. It also verifies that the manifest propagates the source pose,
    ligand template, score, rank, and engine fields.
    """
    parent_workspace = tmp_path / "workspace" / "dock-job"
    current_workspace = tmp_path / "workspace" / "export-job"
    docking_dir = parent_workspace / "docking"
    docking_dir.mkdir(parents=True)

    pose_pdbqt = docking_dir / "dockligand_vina_out.pdbqt"
    ligand_sdf = tmp_path / "inputs" / "ligand.sdf"

    ligand_sdf.parent.mkdir(parents=True)
    pose_pdbqt.write_text("pose pdbqt\n", encoding="utf-8")
    ligand_sdf.write_text("ligand sdf\n", encoding="utf-8")

    parent_manifest_path = docking_dir / "result.json"
    parent_manifest_path.write_text(
        json.dumps(
            {
                "pose_pdbqt": str(pose_pdbqt),
                "ligand_sdf": str(ligand_sdf),
                "score": -7.5,
                "rank": 1,
                "engine": "vina",
            }
        ),
        encoding="utf-8",
    )

    job = FakeJob(current_workspace)
    template_mol = object()
    calls: dict[str, Any] = {}

    monkeypatch.setattr(export_pose, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        export_pose,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )

    def fake_load_first_sdf_molecule(path: Path, *, remove_hs: bool) -> object:
        calls["load_first_sdf_molecule"] = {
            "path": path,
            "remove_hs": remove_hs,
        }
        return template_mol

    def fake_export_pdbqt_to_sdf(
        pdbqt_path: Path,
        output_sdf: Path,
        *,
        template_mol: object,
        add_hydrogens_after_template: bool,
    ) -> Path:
        calls["export_pdbqt_to_sdf"] = {
            "pdbqt_path": pdbqt_path,
            "output_sdf": output_sdf,
            "template_mol": template_mol,
            "add_hydrogens_after_template": add_hydrogens_after_template,
        }
        output_sdf.write_text("exported sdf\n", encoding="utf-8")
        return output_sdf

    monkeypatch.setattr(
        export_pose,
        "load_first_sdf_molecule",
        fake_load_first_sdf_molecule,
    )
    monkeypatch.setattr(
        export_pose,
        "export_pdbqt_to_sdf",
        fake_export_pdbqt_to_sdf,
    )

    export_pose.main(current_workspace)

    export_dir = current_workspace / "export"
    output_sdf = export_dir / "dockligand_vina_out.sdf"
    result_json = export_dir / "result.json"

    assert output_sdf.exists()
    assert result_json.exists()

    assert calls["load_first_sdf_molecule"] == {
        "path": ligand_sdf.resolve(),
        "remove_hs": False,
    }
    assert calls["export_pdbqt_to_sdf"] == {
        "pdbqt_path": pose_pdbqt.resolve(),
        "output_sdf": output_sdf,
        "template_mol": template_mol,
        "add_hydrogens_after_template": True,
    }

    manifest = json.loads(result_json.read_text(encoding="utf-8"))

    assert manifest == {
        "docked_sdf": str(output_sdf.resolve()),
        "engine": "vina",
        "export_dir": str(export_dir),
        "ligand_sdf": str(ligand_sdf.resolve()),
        "pose_pdbqt": str(pose_pdbqt.resolve()),
        "rank": 1,
        "score": -7.5,
    }


def test_export_pose_reports_missing_parent_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail early when the docking-stage manifest is not available.

    A child Grubicy action should not reconstruct upstream paths manually when
    the declared parent product is missing. This test monkeypatches Grubicy's
    ``parent_file`` helper to report that ``docking/result.json`` does not
    exist. The action should fail before loading the ligand template or calling
    the export helper. This keeps missing dependency products visible at the
    workflow boundary.
    """
    job = FakeJob(tmp_path / "workspace" / "export-job")

    monkeypatch.setattr(export_pose, "open_job_from_directory", lambda directory: job)

    def missing_parent_file(
        current_job: Any,
        relative_path: str,
        *,
        must_exist: bool = True,
    ) -> Path:
        raise FileNotFoundError(f"Missing parent file: {relative_path}")

    monkeypatch.setattr(export_pose, "parent_file", missing_parent_file)

    with pytest.raises(FileNotFoundError, match=export_pose.PARENT_MANIFEST):
        export_pose.main(job.workspace)


def test_export_pose_reports_missing_required_parent_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fail when a required file referenced by the parent manifest is missing.

    The export action requires both the docked PDBQT pose and the original
    ligand SDF template. This test provides a parent manifest whose
    ``pose_pdbqt`` points to a missing file while the ligand template exists.
    The action should fail before library conversion starts, giving the user a
    direct path to the missing parent product. Molecular parsing is intentionally
    not involved.
    """
    parent_workspace = tmp_path / "workspace" / "dock-job"
    current_workspace = tmp_path / "workspace" / "export-job"
    docking_dir = parent_workspace / "docking"
    docking_dir.mkdir(parents=True)

    ligand_sdf = tmp_path / "inputs" / "ligand.sdf"
    ligand_sdf.parent.mkdir(parents=True)
    ligand_sdf.write_text("ligand sdf\n", encoding="utf-8")

    parent_manifest_path = docking_dir / "result.json"
    parent_manifest_path.write_text(
        json.dumps(
            {
                "pose_pdbqt": str(docking_dir / "missing_pose.pdbqt"),
                "ligand_sdf": str(ligand_sdf),
            }
        ),
        encoding="utf-8",
    )

    job = FakeJob(current_workspace)

    monkeypatch.setattr(export_pose, "open_job_from_directory", lambda directory: job)
    monkeypatch.setattr(
        export_pose,
        "parent_file",
        lambda current_job, relative_path, *, must_exist=True: parent_manifest_path,
    )

    with pytest.raises(FileNotFoundError, match="pose_pdbqt"):
        export_pose.main(current_workspace)
