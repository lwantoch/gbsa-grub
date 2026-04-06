# gbsa-grub

<p align="center">
  <b>structured molecular workflows</b><br>
  <sub>Docking → MD → GBSA</sub>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/status-dev-orange">
  <img src="https://img.shields.io/badge/version-0.0.1--dev-blue">
  <img src="https://img.shields.io/badge/core-gbsa--pipeline-black">
  <img src="https://img.shields.io/badge/orchestration-GRUBICY-purple">
</p>

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

Phase 1 → connectivity  
(make everything run)

Phase 2 → refinement  
(make it correct)

---

## current status

Docking  ✔  
MD       ✘  
GBSA     ✘  

Goal:

Docking → MD → GBSA (end-to-end)

---

## versioning

0.0.1 → pipeline wired  
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
