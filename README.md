# gbsa-grub

structured molecular workflows  
Docking → MD → GBSA

---

## update (branch: wired)

The project is currently in the *wiring phase*.  
This means that the primary goal is not scientific correctness, but ensuring that all parts of the pipeline are connected and able to exchange data consistently.

At this stage, the system is being shaped into a coherent workflow where each step produces well-defined outputs that can be consumed by the next step. The focus is therefore on structure, reproducibility, and robustness of execution rather than accuracy of results.

---

## what is this?

`gbsa-grub` is an attempt to impose structure on molecular simulation workflows by combining three complementary layers.

The **gbsa-pipeline** provides the scientific building blocks (docking, MD, GBSA).  
**GRUBICY** provides orchestration, i.e. how steps are connected and executed.  
**signac** provides a job-based data model, allowing each computation to be tracked, reproduced, and extended.

Together, they form a system where workflows are no longer implicit scripts, but explicit, inspectable pipelines.

---

## why GRUBICY

In many computational chemistry setups, workflows evolve organically:

run.sh → scripts → files → more scripts → unclear state

This leads to fragile pipelines where:
- dependencies are implicit  
- results are hard to reproduce  
- small changes break downstream steps  

GRUBICY replaces this with a model where:
- each step is a defined action  
- dependencies are explicit  
- execution order is controlled  
- outputs are predictable  

In other words, it turns a collection of scripts into a *system*.

---

## core idea

The pipeline follows a simple conceptual flow:

Docking → MD → GBSA

However, the key idea is not the sequence itself, but the structure behind it.

Each step:
- knows where its inputs come from  
- writes outputs in a defined location  
- can be rerun independently  
- can be extended without breaking the rest  

This enables incremental development without losing control of the workflow.

---

## current status (wired)

| step     | status |
|----------|--------|
| Docking  | ✔ (wired, box currently arbitrary) |
| MD       | ✔ (wired, produces trajectories) |
| GBSA     | ✘ |

At the moment, docking and MD are already connected.  
MD produces trajectories (`.xtc`), which confirms that the pipeline is functionally passing data forward.

However, parameters, physical correctness, and numerical stability are not yet validated.

---

## pipeline view

ligand  
↓  
docking  
↓  
md  
↓  
gbsa  

---

## development strategy

The development follows a strict two-phase philosophy.

**Phase 1 — wiring (current)**  
Everything must run end-to-end.  
Even incorrect outputs are acceptable as long as data flows correctly.

**Phase 2 — refinement**  
Once the pipeline is stable:
- force fields are validated  
- parametrization is corrected  
- physical realism is ensured  

This separation prevents mixing structural bugs with scientific errors.

---

## versioning

| version | meaning |
|--------|--------|
| 0.0.1 | pipeline wired |
| 0.1.0 | docking working |
| 0.2.0 | MD working |
| 0.3.0 | GBSA working |
| 0.4.0 | stabilization |
| 1.0.0 | scientifically usable |

---

## goal

The long-term goal is a fully reproducible workflow:

ligand → docking → MD → GBSA

where each step is:
- transparent  
- reproducible  
- independently testable  

---

## tl;dr

before:
run.sh → ??? → results  

after:
jobs → dependencies → pipeline
