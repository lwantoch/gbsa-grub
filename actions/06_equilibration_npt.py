"""Grubicy runner for unrestrained NPT equilibration.

This action consumes ``equilibration_npt_restrained/result.json`` and runs the
unrestrained NPT equilibration stage. It uses the NPT parameter block without
``-DPOSRES`` so the system can relax after the restrained pressure-equilibration
step. The output is a new ``system.gro``/``system.top`` pair for production.
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
from _md_params import NPT_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_npt_equilibration


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="equilibration_npt")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    equilibrated = run_npt_equilibration(
        100 * BSS.Units.Time.picosecond,
        system,
        work_dir=out_dir,
        params=GromacsParams(**NPT_PARAMS),
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
        params=NPT_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
