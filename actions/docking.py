#!/usr/bin/env python3

from __future__ import annotations

import json
import logging
import shutil
import sys
import traceback
import warnings
from pathlib import Path

import signac
from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    prepare_ligand_with_meeko,
)
from grubicy.typed import load_action_params
from pydantic import BaseModel
from rdkit import Chem

warnings.filterwarnings(
    "ignore",
    message=r".*Redefining '\[electric_potential\]'.*",
)

logging.getLogger("pint.util").setLevel(logging.ERROR)

LOGGER = logging.getLogger(__name__)

if not LOGGER.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )


class DockingParams(BaseModel):
    protein: Path
    ligand: Path


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def ensure_vina() -> None:
    if shutil.which("vina") is None:
        raise RuntimeError("vina binary not found in PATH")


def ensure_obabel() -> None:
    if shutil.which("obabel") is None:
        raise RuntimeError("obabel binary not found in PATH")


def _load_ligand_as_rdkit_mol(ligand_path: Path) -> Chem.Mol:
    ligand_path = Path(ligand_path).resolve()
    suffix = ligand_path.suffix.lower()

    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
        mol = supplier[0] if supplier and len(supplier) > 0 else None
        if mol is None:
            raise ValueError(f"Failed to read ligand SDF: {ligand_path}")
        return mol

    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(ligand_path), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to read ligand MOL: {ligand_path}")
        return mol

    if suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(ligand_path), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to read ligand MOL2: {ligand_path}")
        return mol

    if suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to read ligand PDB: {ligand_path}")
        return mol

    raise ValueError(
        f"Unsupported ligand file type for Meeko preparation: {ligand_path}"
    )


def main(directory: str) -> None:
    project = signac.get_project()
    job = project.open_job(id=Path(directory).name)

    params: DockingParams = load_action_params(job, DockingParams)

    LOGGER.info("Starting docking job %s", job.id)

    ensure_vina()
    ensure_obabel()

    workdir = Path(job.fn("docking"))
    workdir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_docking(
            protein=params.protein.resolve(),
            ligand=params.ligand.resolve(),
            workdir=workdir,
        )

        write_json(workdir / "result.json", result)
        LOGGER.info("Docking finished")

    except Exception as exc:
        LOGGER.error("Docking failed")
        traceback.print_exc()

        write_json(
            workdir / "result.json",
            {
                "status": "failed",
                "error": str(exc),
            },
        )
        raise


def run_docking(protein: Path, ligand: Path, workdir: Path) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)

    engine = VinaEngine(binary="vina", obabel_binary="obabel")

    box = DockingBox(
        center=(7.0, 20.0, 14.0),
        size=(10.0, 10.0, 10.0),
    )

    ligand_pdbqt = workdir / "ligand_input.pdbqt"

    ligand_mol = _load_ligand_as_rdkit_mol(ligand)
    prepare_ligand_with_meeko(
        ligand=ligand_mol,
        output_path=ligand_pdbqt,
        name=ligand.stem,
    )

    request = DockingRequest(
        receptor=protein,
        ligands=[ligand_pdbqt],
        box=box,
        workdir=workdir,
    )

    result = engine.dock(request=request)

    first_pose = result.poses[0] if result.poses else None

    return {
        "status": "success",
        "engine": result.engine,
        "prepared_ligand_pdbqt": str(ligand_pdbqt),
        "final_pose_pdbqt": str(first_pose.pose_path) if first_pose else None,
        "best_score": first_pose.score if first_pose else None,
        "log_file": first_pose.metadata.get("log_file") if first_pose else None,
    }


if __name__ == "__main__":
    main(sys.argv[1])
