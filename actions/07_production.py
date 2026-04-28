"""Grubicy runner for production MD.

This action consumes ``equilibration_npt/result.json`` and runs production with
the production MDP parameter block. The current uploaded ``run_production``
helper does not accept custom parameters, so this runner calls
``run_gro_custom`` directly for production. That keeps the runner compatible
with the uploaded module while still using the same custom-GROMACS mechanism as
the minimization and equilibration helpers.
"""

from __future__ import annotations

import sys

from _common import (
    action_name,
    load_bss_system,
    parent_system_paths,
    save_bss_system,
    stage_dir,
    write_stage_manifest,
)
from _md_params import PRODUCTION_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams, run_gro_custom


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="production")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    production, _protocol = run_gro_custom(
        parameters=GromacsParams(**PRODUCTION_PARAMS),
        system=system,
        work_dir=out_dir,
    )

    output_gro, output_top = save_bss_system(production, out_dir / "system")

    write_stage_manifest(
        job=job,
        stage=stage,
        parent_result=parent_result,
        input_gro=input_gro,
        input_top=input_top,
        output_gro=output_gro,
        output_top=output_top,
        params=PRODUCTION_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
