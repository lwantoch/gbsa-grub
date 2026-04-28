"""Grubicy runner for restrained NVT heating.

This action consumes ``minimization_cg/result.json`` and runs the NVT heating
stage using the custom heating MDP parameter block. The runtime argument is kept
for the public helper signature, but when ``params`` is supplied the current
``run_heating`` helper routes through the custom GROMACS protocol and uses the
runtime encoded by ``dt`` and ``nsteps``. The output is written as a new
``system.gro``/``system.top`` pair for restrained NPT equilibration.
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
from _md_params import HEATING_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_heating


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="equilibration_nvt_heating")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    heated = run_heating(
        50 * BSS.Units.Time.picosecond,
        system,
        work_dir=out_dir,
        params=GromacsParams(**HEATING_PARAMS),
    )

    output_gro, output_top = save_bss_system(heated, out_dir / "system")

    write_stage_manifest(
        job=job,
        stage=stage,
        parent_result=parent_result,
        input_gro=input_gro,
        input_top=input_top,
        output_gro=output_gro,
        output_top=output_top,
        params=HEATING_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
