# /home/grheco/PycharmProjects/gbsa-grub/actions/01_md.py

#!/usr/bin/env python3

from __future__ import annotations

import json
import logging
import sys
import traceback
import warnings
from pathlib import Path
from typing import Any

import signac
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_md_from_docking
from gbsa_pipeline.parametrization import ParametrizationConfig
from gbsa_pipeline.solvation_box import BoxShape, SolvationParams
from pydantic import BaseModel, Field

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


class MDParams(BaseModel):
    """
    Parameters for the gbsa-grub MD action.

    These field names are chosen for workflow readability and Bayesian
    optimization friendliness. They are translated into the concrete
    gbsa_pipeline.solvation_box.SolvationParams field names in
    build_solvation_params().

    Expected statepoint fields
    --------------------------
    protein_pdb:
        Protein PDB to use for parametrization/MD.

    parent_docking_result:
        Relative path inside the parent job that points to the docking result JSON.

    prefer_parent_sdf:
        If True, prefer an SDF from the parent docking result over a PDBQT pose.
        If False, prefer the docked pose first.

    Force-field parameters
    ----------------------
    protein_ff, ligand_ff, charge_method, ligand_net_charge, extra_ff_files

    MD parameters
    -------------
    integrator, nsteps, dt, nstxout_compressed, nstenergy, nstlog

    Solvation parameters
    --------------------
    box_shape:
        Desired solvent box shape ("cubic" or "truncated_octahedron").

    padding_nm:
        Solvent padding shell in nm. This maps to SolvationParams.padding
        and takes precedence over absolute box_size.

    ionic_strength_molar:
        Salt concentration in mol/L. This maps to SolvationParams.ion_concentration.

    neutralize:
        Whether to neutralize the system during solvation.

    Output/control parameters
    -------------------------
    solvate
    save_manifest
    save_parametrized_pdb
    save_solvated_pdb
    save_md_final_pdb
    save_md_native_outputs
    md_output_stem
    """

    protein_pdb: Path
    parent_docking_result: str = "docking/result.json"
    prefer_parent_sdf: bool = True

    protein_ff: str = "ff14SB"
    ligand_ff: str = "gaff2"
    charge_method: str = "am1bcc"
    ligand_net_charge: int | None = None
    extra_ff_files: list[str] = Field(default_factory=list)

    integrator: str = "md"
    nsteps: int = 100
    dt: float = 0.001
    nstxout_compressed: int = 500
    nstenergy: int = 500
    nstlog: int = 500

    box_shape: str = "cubic"
    padding_nm: float = 1.0
    ionic_strength_molar: float = 0.15
    neutralize: bool = True

    solvate: bool = True
    save_manifest: bool = True
    save_parametrized_pdb: bool = True
    save_solvated_pdb: bool = True
    save_md_final_pdb: bool = True
    save_md_native_outputs: bool = True
    md_output_stem: str = "md_final"


def write_json(path: Path, data: dict[str, Any]) -> None:
    """Write JSON data to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    """Read a JSON file from disk."""
    return json.loads(path.read_text(encoding="utf-8"))


def validate_file(path: Path, label: str) -> Path:
    """Validate that a path exists and is a file."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")

    return path


def build_parametrization_config(params: MDParams) -> ParametrizationConfig:
    """Build gbsa-pipeline parametrization config from action params."""
    return ParametrizationConfig(
        protein_ff=params.protein_ff,
        ligand_ff=params.ligand_ff,
        charge_method=params.charge_method,
        extra_ff_files=params.extra_ff_files,
    )


def build_gromacs_params(params: MDParams) -> GromacsParams:
    """Build gbsa-pipeline GROMACS params from action params."""
    return GromacsParams(
        integrator=params.integrator,
        nsteps=params.nsteps,
        dt=params.dt,
        nstxout_compressed=params.nstxout_compressed,
        nstenergy=params.nstenergy,
        nstlog=params.nstlog,
    )


def build_solvation_params(params: MDParams) -> SolvationParams:
    """
    Build gbsa-pipeline solvation params from action params.

    Workflow-facing names are mapped to the actual SolvationParams schema:

    - box_shape -> shape
    - padding_nm -> padding
    - ionic_strength_molar -> ion_concentration
    - neutralize -> neutralize
    """
    return SolvationParams(
        shape=BoxShape(params.box_shape.lower().strip()),
        padding=params.padding_nm,
        ion_concentration=params.ionic_strength_molar,
        neutralize=params.neutralize,
    )


def resolve_pose_from_parent_result(
    result_data: dict[str, Any],
    *,
    prefer_parent_sdf: bool,
) -> tuple[Path, str]:
    """
    Resolve the ligand input for the MD stage from the parent docking result.

    Resolution order
    ----------------
    If prefer_parent_sdf is True:
        1. normalized_ligand_sdf
        2. ligand_input_sdf
        3. final_pose_pdbqt
        4. best_pose_path
        5. docking_result.poses[0].pose_path

    If prefer_parent_sdf is False, the SDF and pose candidate blocks are swapped.

    Returns
    -------
    (resolved_path, source_label)
    """
    sdf_candidates = [
        ("normalized_ligand_sdf", result_data.get("normalized_ligand_sdf")),
        ("ligand_input_sdf", result_data.get("ligand_input_sdf")),
    ]

    pose_candidates = [
        ("final_pose_pdbqt", result_data.get("final_pose_pdbqt")),
        ("best_pose_path", result_data.get("best_pose_path")),
    ]

    docking_result = result_data.get("docking_result", {})
    poses = docking_result.get("poses", [])
    if poses:
        pose_candidates.append(
            ("docking_result.poses[0].pose_path", poses[0].get("pose_path"))
        )

    ordered_candidates = (
        sdf_candidates + pose_candidates
        if prefer_parent_sdf
        else pose_candidates + sdf_candidates
    )

    for source, raw_path in ordered_candidates:
        if not raw_path:
            continue

        candidate = Path(raw_path).expanduser().resolve()
        if candidate.exists() and candidate.is_file():
            return candidate, source

    raise ValueError(
        "Could not resolve ligand input for MD from parent docking result. "
        "Checked SDF/PDBQT candidate fields but found no existing file."
    )


def _native_md_outputs_to_dict(native_outputs: Any) -> dict[str, str | None] | None:
    """
    Convert result.native_md_outputs into a JSON-serializable dictionary.

    This is intentionally defensive because the exact concrete type lives in
    gbsa_pipeline.md and may evolve.
    """
    if native_outputs is None:
        return None

    if hasattr(native_outputs, "to_dict"):
        data = native_outputs.to_dict()
        if isinstance(data, dict):
            return data

    candidate_keys = [
        "run_directory",
        "trajectory_xtc",
        "trajectory_trr",
        "final_gro",
        "portable_run_input_tpr",
        "checkpoint_cpt",
        "energy_edr",
        "md_log",
        "index_ndx",
        "mdp",
    ]

    output: dict[str, str | None] = {}
    found_any = False

    for key in candidate_keys:
        value = getattr(native_outputs, key, None)
        if value is not None:
            found_any = True
            output[key] = str(value)
        else:
            output[key] = None

    return output if found_any or output else None


def main(directory: str) -> None:
    """Run the MD action for one signac job directory."""
    project = signac.get_project()
    job = project.open_job(id=Path(directory).name)

    job_dir = Path(job.path).resolve()
    md_action_dir = Path(job.fn("md_action")).resolve()
    md_action_dir.mkdir(parents=True, exist_ok=True)

    action_result_path = md_action_dir / "result.json"

    try:
        statepoint_path = Path(job.fn("signac_statepoint.json")).resolve()
        statepoint = read_json(statepoint_path)
        params = MDParams.model_validate(statepoint)

        protein_pdb = validate_file(params.protein_pdb, "protein_pdb")

        parent_action = statepoint.get("parent_action")
        if not parent_action:
            raise ValueError(
                "Missing 'parent_action' in signac_statepoint.json for MD job."
            )

        parent_job = project.open_job(id=str(parent_action))
        parent_result_path = validate_file(
            Path(parent_job.fn(params.parent_docking_result)),
            "parent_docking_result",
        )
        parent_result = read_json(parent_result_path)

        docked_ligand_pose, pose_source = resolve_pose_from_parent_result(
            parent_result,
            prefer_parent_sdf=params.prefer_parent_sdf,
        )

        md_work_dir = Path(job.fn("md")).resolve()
        md_work_dir.mkdir(parents=True, exist_ok=True)

        LOGGER.info("====================================================")
        LOGGER.info("MD ACTION")
        LOGGER.info("====================================================")
        LOGGER.info("MD job                 : %s", job.id)
        LOGGER.info("Job dir                : %s", job_dir)
        LOGGER.info("Protein PDB            : %s", protein_pdb)
        LOGGER.info("Parent docking job     : %s", parent_job.id)
        LOGGER.info("Parent result          : %s", parent_result_path)
        LOGGER.info("Selected ligand        : %s", docked_ligand_pose)
        LOGGER.info("Ligand source          : %s", pose_source)
        LOGGER.info("MD work dir            : %s", md_work_dir)
        LOGGER.info("Solvate                : %s", params.solvate)
        LOGGER.info("Save MD final PDB      : %s", params.save_md_final_pdb)
        LOGGER.info("Save native MD outputs : %s", params.save_md_native_outputs)
        LOGGER.info("Solvation shape        : %s", params.box_shape)
        LOGGER.info("Solvation padding (nm) : %s", params.padding_nm)
        LOGGER.info("Ionic strength (M)     : %s", params.ionic_strength_molar)
        LOGGER.info("Neutralize             : %s", params.neutralize)

        result = run_md_from_docking(
            protein_pdb=protein_pdb,
            docked_ligand_pose=docked_ligand_pose,
            work_dir=md_work_dir,
            parametrization_config=build_parametrization_config(params),
            gromacs_params=build_gromacs_params(params),
            solvation_params=build_solvation_params(params),
            ligand_net_charge=params.ligand_net_charge,
            solvate=params.solvate,
            save_manifest=params.save_manifest,
            save_parametrized_pdb=params.save_parametrized_pdb,
            save_solvated_pdb=params.save_solvated_pdb,
            save_md_final_pdb=params.save_md_final_pdb,
            save_md_native_outputs=params.save_md_native_outputs,
            md_output_stem=params.md_output_stem,
        )

        native_md_outputs = _native_md_outputs_to_dict(
            getattr(result, "native_md_outputs", None)
        )

        trajectory_path = None
        portable_tpr_path = None

        if native_md_outputs is not None:
            trajectory_path = native_md_outputs.get(
                "trajectory_xtc"
            ) or native_md_outputs.get("trajectory_trr")
            portable_tpr_path = native_md_outputs.get("portable_run_input_tpr")

        write_json(
            action_result_path,
            {
                "status": "success",
                "job_id": job.id,
                "job_dir": str(job_dir),
                "protein_pdb": str(protein_pdb),
                "parent_job_id": parent_job.id,
                "parent_result_path": str(parent_result_path),
                "selected_pose_source": pose_source,
                "selected_pose_path": str(docked_ligand_pose),
                "md_work_dir": str(md_work_dir),
                "md_manifest_path": str(result.manifest_path)
                if result.manifest_path
                else None,
                "prepared_protein_pdb": str(result.prepared_protein_pdb),
                "normalized_ligand_sdf": str(result.normalized_ligand_sdf),
                "parametrization_dir": str(result.parametrization_dir),
                "solvation_dir": str(result.solvation_dir),
                "md_dir": str(result.md_dir),
                "parametrized_pdb": str(result.parametrized_pdb)
                if result.parametrized_pdb is not None
                else None,
                "solvated_pdb": str(result.solvated_pdb)
                if result.solvated_pdb is not None
                else None,
                "final_md_pdb": str(result.final_md_pdb)
                if result.final_md_pdb is not None
                else None,
                "trajectory_path": trajectory_path,
                "tpr_path": portable_tpr_path,
                "native_md_outputs": native_md_outputs,
            },
        )

        LOGGER.info("MD action finished successfully")
        if trajectory_path is not None:
            LOGGER.info("Trajectory             : %s", trajectory_path)
        if portable_tpr_path is not None:
            LOGGER.info("TPR                    : %s", portable_tpr_path)

    except Exception as exc:
        LOGGER.error("MD action failed")
        traceback.print_exc()

        write_json(
            action_result_path,
            {
                "status": "failed",
                "job_id": job.id,
                "job_dir": str(job_dir),
                "error": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            },
        )
        raise


if __name__ == "__main__":
    main(sys.argv[1])
