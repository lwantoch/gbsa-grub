"""Grubicy runner for steepest-descent minimization.

This action consumes ``system_setup/result.json`` and runs the first GROMACS
energy minimization stage using the SD parameter block copied from the
integration-test workflow. It writes a new disk handoff system into
``minimization_sd/system.gro`` and ``minimization_sd/system.top``. The next
action should consume only this action's ``result.json``.
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
from _md_params import SD_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_minimization


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="minimization_sd")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    minimized = run_minimization(
        system,
        work_dir=out_dir,
        params=GromacsParams(**SD_PARAMS),
    )

    output_gro, output_top = save_bss_system(minimized, out_dir / "system")

    write_stage_manifest(
        job=job,
        stage=stage,
        parent_result=parent_result,
        input_gro=input_gro,
        input_top=input_top,
        output_gro=output_gro,
        output_top=output_top,
        params=SD_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
