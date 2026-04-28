"""Grubicy runner for parametrization and solvation.

This action consumes ``docking/result.json`` and creates the first MD handoff
system. It parameterizes the protein and docked ligand, then immediately
solvates the parametrized complex with OpenMM/ParmEd. Parametrization and
solvation stay in one action because the solvation helper currently reuses the
in-memory force field and ParmEd structure from parametrization. The action
writes ``system_setup/system.gro``, ``system_setup/system.top``, and
``system_setup/result.json`` for the first minimization action.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from _common import action_name, parent_manifest, stage_dir, write_json
from gbsa_pipeline.parametrization import ParametrizationInput, parametrize
from gbsa_pipeline.solvation_box import SolvationParams
from gbsa_pipeline.solvation_openmm import solvate_openmm


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)

    stage = action_name(job, fallback="system_setup")
    out_dir = stage_dir(job, stage)

    _parent, parent_result_path, docking_manifest = parent_manifest(job)

    protein_pdb = Path(str(docking_manifest["protein_pdb"]))
    docked_ligand_sdf = Path(str(docking_manifest["docked_ligand_sdf"]))

    parametrization_dir = out_dir / "parametrization"
    solvation_dir = out_dir / "solvation"
    parametrization_dir.mkdir(parents=True, exist_ok=True)
    solvation_dir.mkdir(parents=True, exist_ok=True)

    parametrized = parametrize(
        ParametrizationInput(
            protein_pdb=protein_pdb,
            ligand_sdf=docked_ligand_sdf,
            work_dir=parametrization_dir,
        )
    )

    solvated = solvate_openmm(
        parametrized=parametrized,
        params=SolvationParams(
            water_model="tip3p",
            shape="truncated_octahedron",
            padding=1.0,
            ion_concentration=0.15,
            neutralize=True,
        ),
        output_gro=solvation_dir / "solvated.gro",
        output_top=solvation_dir / "solvated.top",
    )

    output_gro = out_dir / "system.gro"
    output_top = out_dir / "system.top"

    shutil.copyfile(solvated.gro_file, output_gro)
    shutil.copyfile(solvated.top_file, output_top)

    write_json(
        out_dir / "result.json",
        {
            "stage": stage,
            "parent_result": parent_result_path,
            "protein_pdb": protein_pdb,
            "docked_ligand_sdf": docked_ligand_sdf,
            "parametrized_gro": parametrized.gro_file,
            "parametrized_top": parametrized.top_file,
            "crystal_waters_pdb": parametrized.crystal_waters_pdb,
            "solvated_gro": solvated.gro_file,
            "solvated_top": solvated.top_file,
            "output_gro": output_gro,
            "output_top": output_top,
            "system_gro": output_gro,
            "system_top": output_top,
        },
    )


if __name__ == "__main__":
    main(sys.argv[1])
