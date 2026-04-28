"""Grubicy runner for docking.

This runner follows the docking part of the inspectable integration test. It
loads the ligand SDF, prepares ligand PDBQT with Meeko, prepares receptor PDBQT,
runs Vina, exports the docked pose back to SDF with template bond-order repair,
and writes ``docking/result.json``. The next action should consume only that
manifest, not guess Vina output names. The required state point keys are
``protein_pdb``, ``ligand_sdf``, ``box_center``, and ``box_size``.
"""

from __future__ import annotations

import sys
from typing import Any

from _common import action_name, require_sp_path, stage_dir, statepoint, write_json
from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    convert_receptor_pdb_to_pdbqt,
    export_pdbqt_to_sdf,
    load_first_sdf_molecule,
    prepare_ligand_with_meeko,
)

MEEKO_RECEPTOR_BINARY = "mk_prepare_receptor.py"


def _float_tuple(value: Any, *, length: int, name: str) -> tuple[float, ...]:
    if not isinstance(value, (list, tuple)):
        raise TypeError(f"{name} must be a list or tuple.")
    if len(value) != length:
        raise ValueError(f"{name} must contain exactly {length} values.")
    return tuple(float(item) for item in value)


def _vina_parameters(sp: dict[str, Any]) -> dict[str, Any]:
    if isinstance(sp.get("vina_parameters"), dict):
        return dict(sp["vina_parameters"])

    allowed = {"exhaustiveness", "num_modes", "energy_range", "cpu", "verbosity"}
    return {key: sp[key] for key in allowed if key in sp}


def main(directory: str) -> None:
    from _common import (
        open_job_from_directory,  # local import keeps script execution explicit
    )

    job = open_job_from_directory(directory)
    sp = statepoint(job)

    stage = action_name(job, fallback="docking")
    out_dir = stage_dir(job, stage)

    protein_pdb = require_sp_path(job, "protein_pdb")
    ligand_sdf = require_sp_path(job, "ligand_sdf")

    box = DockingBox(
        center=_float_tuple(sp["box_center"], length=3, name="box_center"),
        size=_float_tuple(sp["box_size"], length=3, name="box_size"),
    )

    ligand_name = str(sp.get("ligand_id", ligand_sdf.stem))
    ligand_pdbqt = out_dir / f"{ligand_name}.pdbqt"
    receptor_pdbqt = out_dir / "receptor.pdbqt"
    docked_sdf = out_dir / "docked_ligand.sdf"

    ligand_molecule = load_first_sdf_molecule(ligand_sdf, remove_hs=False)

    prepared_ligand = prepare_ligand_with_meeko(
        ligand_molecule,
        ligand_pdbqt,
        name=ligand_name,
    )

    prepared_receptor = convert_receptor_pdb_to_pdbqt(
        protein_pdb,
        output_path=receptor_pdbqt,
        mk_prepare_receptor_binary=str(
            sp.get("mk_prepare_receptor_binary", MEEKO_RECEPTOR_BINARY)
        ),
    )

    engine = VinaEngine(
        binary=str(sp.get("vina_binary", "vina")),
        mk_prepare_receptor_binary=str(
            sp.get("mk_prepare_receptor_binary", MEEKO_RECEPTOR_BINARY)
        ),
    )

    request = DockingRequest(
        receptor=prepared_receptor,
        ligands=[prepared_ligand],
        box=box,
        seed=sp.get("seed"),
        workdir=out_dir,
        parameters=_vina_parameters(sp),
    )

    docking_result = engine.dock(request=request)
    if not docking_result.poses:
        raise RuntimeError("Vina returned no docked poses.")

    best_pose = docking_result.poses[0]

    export_pdbqt_to_sdf(
        best_pose.pose_path,
        docked_sdf,
        template_mol=ligand_molecule,
        add_hydrogens_after_template=True,
    )

    write_json(
        out_dir / "result.json",
        {
            "stage": stage,
            "protein_pdb": protein_pdb,
            "ligand_template_sdf": ligand_sdf,
            "prepared_ligand_pdbqt": prepared_ligand,
            "prepared_receptor_pdbqt": prepared_receptor,
            "docked_pose_pdbqt": best_pose.pose_path,
            "docked_ligand_sdf": docked_sdf,
            "score": best_pose.score,
            "rank": best_pose.rank,
            "engine": docking_result.engine,
            "parameters": dict(docking_result.parameters),
            "pose_metadata": dict(best_pose.metadata),
        },
    )


if __name__ == "__main__":
    main(sys.argv[1])
