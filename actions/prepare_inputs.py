"""Prepare receptor and ligand inputs for the grubicy docking workflow.

This action is the root step of the first gbsa-grub smoke workflow. It reads the
protein PDB path, ligand SDF path, and docking box values from the Grubicy
action parameters stored in the signac statepoint. It then delegates all
chemistry-specific work to ``gbsa_pipeline.docking`` by loading the ligand
molecule, preparing ligand PDBQT with Meeko, and preparing receptor PDBQT
through the library helper. The action writes only declared workflow products
under ``prepare/`` so downstream Grubicy actions can consume stable files rather
than in-memory Python objects.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Any

from gbsa_pipeline.docking import (
    DockingBox,
    convert_receptor_pdb_to_pdbqt,
    load_first_sdf_molecule,
    prepare_ligand_with_meeko,
)
from grubicy import open_job_from_directory
from grubicy.typed import load_action_params
from pydantic import BaseModel, ConfigDict

MEEKO_RECEPTOR_BINARY = "mk_prepare_receptor.py"


class PrepareInputsParams(BaseModel):
    """Validated parameters for the prepare action.

    Grubicy materializes the action parameters into the signac statepoint, and
    ``load_action_params`` validates those values against this model at runtime.
    The box fields are fixed-length three-float tuples so malformed TOML arrays
    fail before any chemistry helper is called. Paths are validated here only as
    path-like values; existence and file-kind checks remain explicit in
    ``_required_input_file`` so the error messages point to the workflow input
    being inspected. Extra statepoint keys are forbidden because this action
    should consume only the parameters declared for the prepare step.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    protein_pdb: Path
    ligand_sdf: Path
    box_center: tuple[float, float, float]
    box_size: tuple[float, float, float]


def _required_input_file(path: Path, *, key: str) -> Path:
    """Return a required input path as a resolved file path.

    The prepare action is the root file-ingestion step for the smoke workflow,
    so missing input paths should fail before any chemistry helper is called.
    This helper keeps file-existence and file-kind checks at the workflow
    boundary. It does not inspect molecular contents; ligand and receptor
    parsing remain delegated to ``gbsa_pipeline.docking``. The resolved absolute
    path is written into the prepare-stage manifest.
    """
    resolved = path.expanduser().resolve()

    if not resolved.exists():
        raise FileNotFoundError(
            f"Input file not found for parameter '{key}': {resolved}"
        )

    if not resolved.is_file():
        raise ValueError(f"Input path for parameter '{key}' is not a file: {resolved}")

    return resolved


def _write_json(path: Path, data: dict[str, Any]) -> Path:
    """Write a small JSON manifest with deterministic formatting.

    The prepare action writes ``prepare/result.json`` as the contract consumed by
    the docking action. The manifest records the prepared ligand and receptor
    paths plus the original input paths and box values. The JSON file should
    contain only lightweight metadata and file paths; molecular structures are
    written as normal workflow products. Serialization errors are not hidden
    because non-JSON-compatible manifest values are runner bugs.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def main(directory: str | Path) -> None:
    """Prepare ligand PDBQT, receptor PDBQT, and a prepare-stage manifest.

    The action opens the current job from the workspace directory supplied by
    row/Grubicy. Action parameters are loaded through Grubicy's typed runtime
    helper instead of being pulled manually from ``job.sp``. Chemistry-specific
    work remains in ``gbsa_pipeline.docking`` so the action runner only
    coordinates input files, output paths, and the declared prepare manifest.
    The resulting files are written under ``prepare/`` with stable names for
    downstream actions.
    """
    job = open_job_from_directory(str(directory))
    params = load_action_params(job, PrepareInputsParams)

    if shutil.which(MEEKO_RECEPTOR_BINARY) is None:
        raise RuntimeError(f"{MEEKO_RECEPTOR_BINARY} not available in PATH.")

    protein_pdb = _required_input_file(params.protein_pdb, key="protein_pdb")
    ligand_sdf = _required_input_file(params.ligand_sdf, key="ligand_sdf")

    box = DockingBox(
        center=params.box_center,
        size=params.box_size,
    )

    prepare_dir = Path(job.fn("prepare"))
    prepare_dir.mkdir(parents=True, exist_ok=True)

    ligand_pdbqt = prepare_dir / "ligand.pdbqt"
    receptor_pdbqt = prepare_dir / "receptor.pdbqt"
    result_json = prepare_dir / "result.json"

    ligand_molecule = load_first_sdf_molecule(ligand_sdf, remove_hs=False)

    prepared_ligand = prepare_ligand_with_meeko(
        ligand_molecule,
        ligand_pdbqt,
        name="DOCKLIG",
    )

    prepared_receptor = convert_receptor_pdb_to_pdbqt(
        protein_pdb,
        output_path=receptor_pdbqt,
        mk_prepare_receptor_binary=MEEKO_RECEPTOR_BINARY,
    )

    _write_json(
        result_json,
        {
            "protein_pdb": str(protein_pdb),
            "ligand_sdf": str(ligand_sdf),
            "ligand_pdbqt": str(Path(prepared_ligand).resolve()),
            "receptor_pdbqt": str(Path(prepared_receptor).resolve()),
            "box_center": list(box.center),
            "box_size": list(box.size),
            "prepare_dir": str(prepare_dir),
            "mk_prepare_receptor_binary": MEEKO_RECEPTOR_BINARY,
        },
    )


if __name__ == "__main__":
    main(sys.argv[1])
