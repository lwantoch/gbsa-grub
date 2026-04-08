# gbsa-grub

<p align="center">
  <b>structured molecular workflows</b><br>
  <sub>Docking → MD → GBSA</sub>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/status-dev-orange">
  <img src="https://img.shields.io/badge/version-0.0.1--wired-blue">
  <img src="https://img.shields.io/badge/core-gbsa--pipeline-black">
  <img src="https://img.shields.io/badge/orchestration-GRUBICY-purple">
</p>

---

## update (branch: wired)

Current branch: **wired**

Status:
- full pipeline connectivity is being established
- all modules are being connected end-to-end
- correctness is NOT yet the priority
- structure and data flow are the focus

Key changes:

- GRUBICY pipeline (`pipeline.toml`) defined and validated
- action naming unified (`00_docking`, `01_md`, ...)
- dependency chain established via `deps`
- signac-compatible job structure in place
- MD stage connected to docking outputs
- unified manifest/result.json handling across steps
- `unigbsa.py` refactored to be workflow-compatible:
  - request/result model introduced
  - deterministic output structure
  - JSON serialization for pipeline integration

Important principle:

> First make everything run. Then make it correct.

---

## what is this?

`gbsa-grub` connects:

- gbsa-pipeline → scientific logic  
- GRUBICY → workflow structure  
- signac → job model  

Goal: structured, reproducible molecular workflows.

---

## why GRUBICY

Without:

scripts → files → chaos

With:

jobs → dependencies → pipeline

GRUBICY gives:

- explicit parent-child jobs  
- reproducible pipelines  
- central pipeline definition  
- safe schema evolution  

---

## core idea

Docking → MD → GBSA

Each step:

- knows its parent  
- stores its own data  
- is reproducible  

---

## development strategy

Phase 1 → connectivity (current)
- wire everything
- ensure data flows through all stages
- tolerate incorrect science

Phase 2 → refinement
- fix chemistry
- fix parametrization
- validate outputs

---

## current status (wired)

Docking  ✔ (connected, produces outputs)
MD       ✔ (connected, not yet validated)
GBSA     ✘ (next step)

Goal:

Docking → MD → GBSA (end-to-end execution)

---

## versioning

0.0.1 → pipeline wired (current)  
0.1.0 → docking working  
0.2.0 → MD working  
0.3.0 → GBSA working  
0.4.0 → stabilization  
1.0.0 → scientifically usable  

---

## goal

ligand → docking → MD → GBSA

---

## tl;dr

before:
run.sh → ??? → results

after:
jobs → dependencies → pipeline
