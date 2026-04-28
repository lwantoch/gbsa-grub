# /home/grheco/repositorios/gbsa-pipeline/actions/03_minimization_cg.py

"""Grubicy runner for conjugate-gradient minimization.

This action consumes ``minimization_sd/result.json`` and runs the second
minimization step with the CG parameter block from the integration-test
workflow. It does not repeat parametrization, solvation, or docking. Its only
job is to load the parent ``system.gro``/``system.top``, run CG minimization,
and write a new handoff system for the NVT heating action.
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
from _md_params import CG_PARAMS
from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.md import run_minimization


def main(directory: str) -> None:
    from _common import open_job_from_directory

    job = open_job_from_directory(directory)
    stage = action_name(job, fallback="minimization_cg")
    out_dir = stage_dir(job, stage)

    input_gro, input_top, parent_result, _manifest = parent_system_paths(job)
    system = load_bss_system(input_gro, input_top)

    minimized = run_minimization(
        system,
        work_dir=out_dir,
        params=GromacsParams(**CG_PARAMS),
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
        params=CG_PARAMS,
    )


if __name__ == "__main__":
    main(sys.argv[1])
