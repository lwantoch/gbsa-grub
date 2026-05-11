# gbsa-grub

Structured molecular workflows  
Docking → MD → GBSA

---

## status

![status](https://img.shields.io/badge/status-workflow--skeleton-blue)
![version](https://img.shields.io/badge/version-0.0.1-lightgrey)
![phase](https://img.shields.io/badge/phase-orchestration-important)
![orchestration](https://img.shields.io/badge/orchestration-GRUBICY-black)
![python](https://img.shields.io/badge/python-3.12+-blue)
![dev](https://img.shields.io/badge/state-under--development-orange)

**Current state:** first workflow skeleton for Docking → MD → GBSA development.

This repository defines the orchestration layer for a molecular simulation workflow built around `gbsa-pipeline`. The current focus is to make the workflow explicit: each stage has declared inputs, declared outputs, parent dependencies, and job-local artifacts that can be inspected independently.

---

## what is this?

`gbsa-grub` connects the computational workflow

```text
ligand/protein input → docking → MD → GBSA
```

to a structured action-based execution model.

The repository combines three layers:

| layer | role |
|------|------|
| **gbsa-pipeline** | scientific library modules: docking, parametrization, solvation, MD, GBSA |
| **GRUBICY** | workflow orchestration: actions, dependencies, declared products |
| **signac** | job and workspace model |

`gbsa-grub` itself should remain thin. It should not reimplement chemistry or simulation logic. Its job is to call `gbsa-pipeline` modules in the correct order and make every workflow boundary visible on disk.

---

## workflow overview

Current action chain:

```text
prepare_inputs
→ dock
→ export_pose
→ setup_system
→ md_sd
→ md_cg
→ md_nvt_res
→ md_npt_res
→ md_npt
→ md_production
→ gbsa
```

The current implementation wires the Docking → MD path and prepares the structure for adding GBSA as the next workflow action. Longer MD calculations and the GBSA stage are intended to run on CESDGA once the MD products are generated in the expected environment.

---

## actions

### `prepare_inputs`

Prepares docking inputs from the initial protein and ligand files.

Declared products:

```text
prepare/ligand.pdbqt
prepare/receptor.pdbqt
prepare/result.json
```

### `dock`

Runs AutoDock Vina using the prepared ligand and receptor.

Declared products:

```text
docking/dockligand_vina_out.pdbqt
docking/dockligand_vina.log
docking/result.json
```

The docking PDBQT is treated as a pose container. It is not the final ligand chemistry source for parametrization.

### `export_pose`

Exports the docked ligand pose to SDF for downstream parametrization.

Declared products:

```text
export/dockligand_vina_out.sdf
export/result.json
```

The exported SDF is the file consumed by the system setup step.

### `setup_system`

Runs parametrization and solvation and writes the solvated GROMACS system.

Declared products:

```text
setup/solvated.gro
setup/solvated.top
setup/result.json
```

This action currently keeps parametrization and solvation together because the solvation helper uses in-memory objects returned by parametrization.

### MD stages

All MD actions use the shared runner:

```text
actions/md_stage.py
```

Each MD stage reads the parent `.gro/.top`, runs one MD helper from `gbsa-pipeline`, saves a new `.gro/.top`, and writes a `result.json`.

Declared products:

```text
sd/system.gro
sd/system.top
sd/result.json

cg/system.gro
cg/system.top
cg/result.json

nvt_res/system.gro
nvt_res/system.top
nvt_res/result.json

npt_res/system.gro
npt_res/system.top
npt_res/result.json

npt/system.gro
npt/system.top
npt/result.json

production/system.gro
production/system.top
production/result.json
```

### `gbsa`

GBSA is the next workflow target after MD production outputs are available.

Planned role:

```text
production MD outputs → GBSA input preparation → endpoint scoring
```

The intended GBSA implementation belongs to `gbsa-pipeline`; `gbsa-grub` will provide the action boundary and file handoff.

---

## design

The repository follows a thin-runner design.

`gbsa-grub` owns:

```text
- pipeline.toml
- grubicy action definitions
- action runners
- parent/child file handoffs
- job-local result manifests
```

`gbsa-pipeline` owns:

```text
- docking helpers
- ligand/receptor preparation
- PDBQT → SDF export and chemistry reconstruction
- parametrization
- solvation
- MD helpers
- GBSA functionality
```

This separation keeps the scientific implementation reusable outside grubicy and keeps the workflow layer easy to inspect.

---

## repository layout

```text
gbsa-grub/
├── actions/
│   ├── prepare_inputs.py
│   ├── dock.py
│   ├── export_pose.py
│   ├── setup_system.py
│   └── md_stage.py
├── config/
│   └── md_params.py
├── data/
│   ├── protein.pdb
│   └── ligand.sdf
├── tests/
│   └── unit/
├── pipeline.toml
└── workspace/
```

`workspace/` contains generated job directories and should not be committed.

---

## workflow specification

The workflow is defined in:

```text
pipeline.toml
```

The first smoke experiment currently uses:

```text
data/protein.pdb
data/ligand.sdf
```

Each action declares its products explicitly so grubicy can track workflow progress through file existence and parent dependencies.

---

## testing

Static workflow tests check the `pipeline.toml` contract:

```bash
pixi run pytest tests/unit/test_pipeline_toml.py -q
```

Action-level unit tests check the thin runners without running expensive external tools:

```bash
pixi run pytest tests/unit/test_prepare_inputs_action.py -q
pixi run pytest tests/unit/test_dock_action.py -q
pixi run pytest tests/unit/test_export_pose_action.py -q
pixi run pytest tests/unit/test_setup_system_action.py -q
pixi run pytest tests/unit/test_md_stage_action.py -q
```

Validate the grubicy configuration:

```bash
pixi run grubicy validate pipeline.toml
```

Inspect the next runnable action:

```bash
pixi run grubicy submit pipeline.toml --dry-run
```

Submit one runnable action:

```bash
pixi run grubicy submit pipeline.toml --limit 1
```

---

## development phases

| version | target |
|--------|--------|
| 0.0.1 | workflow skeleton wired |
| 0.1.0 | docking stage stabilized |
| 0.2.0 | MD stage stabilized |
| 0.3.0 | GBSA stage added and working |
| 0.4.0 | workflow stabilization |
| 1.0.0 | scientifically usable workflow |

---

## execution environment

Short smoke tests can be run locally.

Longer MD and GBSA-oriented runs are expected to move to CESDGA because production trajectories and intermediate MD artifacts can exceed local disk capacity.

---

## tl;dr

`gbsa-grub` is the orchestration layer for a Docking → MD → GBSA workflow.

```text
gbsa-pipeline = scientific library
gbsa-grub     = workflow structure
GRUBICY       = action orchestration
signac        = job/workspace model
```

The goal is a workflow where each stage is explicit, inspectable, and reusable as part of a larger molecular simulation pipeline.
