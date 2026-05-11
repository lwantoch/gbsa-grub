"""Static tests for the grubicy pipeline specification.

These tests validate the workflow specification before any grubicy command is
started. They intentionally do not run grubicy, row, signac, docking, or MD,
because the goal is to catch simple wiring mistakes in ``pipeline.toml`` early.
The tested properties are the stable orchestration contract: action order,
parent dependencies, declared output products, runner commands, and the first
smoke-test experiment. Runtime behavior is still covered by running
``grubicy validate`` and ``grubicy submit --dry-run`` outside this unit test.
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

PIPELINE_TOML = Path(__file__).resolve().parents[2] / "pipeline.toml"

EXPECTED_ACTIONS = [
    "prepare_inputs",
    "dock",
    "export_pose",
    "setup_system",
    "md_sd",
    "md_cg",
    "md_nvt_res",
    "md_npt_res",
    "md_npt",
    "md_production",
]

EXPECTED_DEPS = {
    "dock": "prepare_inputs",
    "export_pose": "dock",
    "setup_system": "export_pose",
    "md_sd": "setup_system",
    "md_cg": "md_sd",
    "md_nvt_res": "md_cg",
    "md_npt_res": "md_nvt_res",
    "md_npt": "md_npt_res",
    "md_production": "md_npt",
}

EXPECTED_OUTPUTS = {
    "prepare_inputs": [
        "prepare/ligand.pdbqt",
        "prepare/receptor.pdbqt",
        "prepare/result.json",
    ],
    "dock": [
        "docking/dockligand_vina_out.pdbqt",
        "docking/dockligand_vina.log",
        "docking/result.json",
    ],
    "export_pose": [
        "export/dockligand_vina_out.sdf",
        "export/result.json",
    ],
    "setup_system": [
        "setup/solvated.gro",
        "setup/solvated.top",
        "setup/result.json",
    ],
    "md_sd": [
        "sd/system.gro",
        "sd/system.top",
        "sd/result.json",
    ],
    "md_cg": [
        "cg/system.gro",
        "cg/system.top",
        "cg/result.json",
    ],
    "md_nvt_res": [
        "nvt_res/system.gro",
        "nvt_res/system.top",
        "nvt_res/result.json",
    ],
    "md_npt_res": [
        "npt_res/system.gro",
        "npt_res/system.top",
        "npt_res/result.json",
    ],
    "md_npt": [
        "npt/system.gro",
        "npt/system.top",
        "npt/result.json",
    ],
    "md_production": [
        "production/system.gro",
        "production/system.top",
        "production/result.json",
    ],
}

EXPECTED_RUNNERS = {
    "prepare_inputs": "pixi run python actions/prepare_inputs.py {directory}",
    "dock": "pixi run python actions/dock.py {directory}",
    "export_pose": "pixi run python actions/export_pose.py {directory}",
    "setup_system": "pixi run python actions/setup_system.py {directory}",
    "md_sd": "pixi run python actions/md_stage.py {directory}",
    "md_cg": "pixi run python actions/md_stage.py {directory}",
    "md_nvt_res": "pixi run python actions/md_stage.py {directory}",
    "md_npt_res": "pixi run python actions/md_stage.py {directory}",
    "md_npt": "pixi run python actions/md_stage.py {directory}",
    "md_production": "pixi run python actions/md_stage.py {directory}",
}


def _load_pipeline() -> dict[str, Any]:
    """Load the grubicy TOML specification as a plain dictionary.

    The test uses Python's standard ``tomllib`` parser so the file is checked as
    TOML before grubicy receives it. This keeps the unit test independent from
    grubicy internals and avoids materializing real workflow jobs. The returned
    dictionary is intentionally not wrapped in custom models because the test is
    only checking a small stable contract. Deeper semantic validation still
    belongs to ``grubicy validate``.
    """
    with PIPELINE_TOML.open("rb") as handle:
        return tomllib.load(handle)


def test_pipeline_toml_contains_expected_action_chain() -> None:
    """Check that the workflow action chain is present in the expected order.

    The order matters because this file describes a linear smoke workflow from
    input preparation to production MD. Accidentally renaming or reordering an
    action would break parent dependency assumptions in the runner scripts. The
    test also checks uniqueness so duplicate action names are caught before
    grubicy materialization. This is a static contract test and does not imply
    that the action scripts themselves run successfully.
    """
    pipeline = _load_pipeline()
    action_names = [action["name"] for action in pipeline["actions"]]

    assert action_names == EXPECTED_ACTIONS
    assert len(action_names) == len(set(action_names))


def test_pipeline_toml_contains_expected_parent_dependencies() -> None:
    """Check that every child action depends on the expected parent action.

    Grubicy uses declared dependencies to connect child jobs to parent products.
    The first action should not have a parent, while every later action should
    point to exactly the previous workflow stage through ``parent_action``. This
    test catches broken or misspelled dependency names before runtime. It does
    not inspect files in the workspace because no jobs are materialized here.
    """
    pipeline = _load_pipeline()
    actions = {action["name"]: action for action in pipeline["actions"]}

    assert "deps" not in actions["prepare_inputs"]

    for action_name, parent_action in EXPECTED_DEPS.items():
        deps = actions[action_name]["deps"]

        assert deps["action"] == parent_action
        assert deps["sp_key"] == "parent_action"


def test_pipeline_toml_declares_expected_outputs() -> None:
    """Check that each action declares the expected workflow products.

    Grubicy tracks completion through declared output files, so these paths are
    part of the workflow contract. Every action should write a small
    ``result.json`` manifest plus the concrete scientific artefacts needed by
    downstream steps. Testing exact paths catches naming drift between the TOML
    spec and the action scripts. The files do not need to exist for this static
    test.
    """
    pipeline = _load_pipeline()
    actions = {action["name"]: action for action in pipeline["actions"]}

    for action_name, outputs in EXPECTED_OUTPUTS.items():
        assert actions[action_name]["outputs"] == outputs


def test_pipeline_toml_declares_expected_runners() -> None:
    """Check that each action uses the expected runner command.

    The workflow intentionally uses one shared ``actions/md_stage.py`` runner for
    all MD stages while the earlier stages each have their own thin action file.
    This test catches accidental runner drift, especially for the MD stages where
    duplicate runner files should not be introduced. The command strings are
    tested statically because the actual execution belongs to grubicy dry-run or
    integration checks. Keeping this contract explicit makes the TOML easier to
    review.
    """
    pipeline = _load_pipeline()
    actions = {action["name"]: action for action in pipeline["actions"]}

    for action_name, runner in EXPECTED_RUNNERS.items():
        assert actions[action_name]["runner"] == runner


def test_pipeline_toml_has_first_smoke_experiment() -> None:
    """Check that the first smoke-test experiment provides root parameters.

    The first experiment should provide the local input files and docking box for
    ``prepare_inputs`` plus the solvent settings for ``setup_system``. Child MD
    actions do not need their own experiment parameters yet because their
    current parameter blocks live in the runner layer. This test keeps the
    initial example complete enough for ``grubicy prepare`` to materialize the
    workflow. It intentionally does not check that the referenced files exist,
    because this static test should not depend on local data checkout state.
    """
    pipeline = _load_pipeline()
    experiments = pipeline["experiment"]

    assert len(experiments) == 1

    experiment = experiments[0]

    assert experiment["prepare_inputs"]["protein_pdb"] == "data/protein.pdb"
    assert experiment["prepare_inputs"]["ligand_sdf"] == "data/ligand.sdf"
    assert experiment["prepare_inputs"]["box_center"] == [10.115, 39.148, 53.112]
    assert experiment["prepare_inputs"]["box_size"] == [10.0, 10.0, 10.0]

    assert experiment["setup_system"]["water_model"] == "tip3p"
    assert experiment["setup_system"]["box_shape"] == "truncated_octahedron"
    assert experiment["setup_system"]["padding_nm"] == 1.0
    assert experiment["setup_system"]["ion_concentration"] == 0.15
    assert experiment["setup_system"]["neutralize"] is True
