"""Grubicy runner for restrained NPT equilibration.

This action consumes ``equilibration_nvt_heating/result.json`` and runs the
restrained NPT equilibration stage. The restraint behavior comes from the custom
MDP block through ``define = -DPOSRES`` and ``refcoord_scaling = com``. The
runner does not know chemistry or topology details; it only runs one MD stage
from the parent disk handoff system and writes the next disk handoff system.
"""

from __future__ import annotations

import sys

import BioSimSpace as BSS
from _common import (
    action_name,
    load_bss_system,
    parent_system_paths,
    save_bss_system,
    stage_dir,
    write_stage_manifest,
)
from _md_params import NPT_RESTRAINED_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_npt_equilibration


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="equilibration_npt_restrained")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    equilibrated = run_npt_equilibration(
        100 * BSS.Units.Time.picosecond,
        system,
        work_dir=out_dir,
        params=GromacsParams(**NPT_RESTRAINED_PARAMS),
    )

    output_gro, output_top = save_bss_system(equilibrated, out_dir / "system")

    write_stage_manifest(
        job=job,
        stage=stage,
        parent_result=parent_result,
        input_gro=input_gro,
        input_top=input_top,
        output_gro=output_gro,
        output_top=output_top,
        params=NPT_RESTRAINED_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
