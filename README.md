# gbsa-grub

Structured molecular workflows  
Docking → MD → GBSA

---

## status (wired branch)

![status](https://img.shields.io/badge/status-wired-blue)
![version](https://img.shields.io/badge/version-0.0.1-lightgrey)
![phase](https://img.shields.io/badge/phase-wiring-important)
![orchestration](https://img.shields.io/badge/orchestration-GRUBICY-black)
![python](https://img.shields.io/badge/python-3.11+-blue)
![dev](https://img.shields.io/badge/state-under--development-orange)

**Active branch:** `wired`  
**Current state:** first end-to-end version connected, not yet ready for routine scientific use

This repository currently lives in its **wiring phase**. That means the central question is not yet *“Are the numbers trustworthy?”* but rather *“Can the workflow hold its own shape?”* At this stage, the project is about forcing a molecular simulation pipeline into something explicit: well-defined inputs, well-defined outputs, clear dependencies, and a job structure that can be inspected instead of guessed.

---

## what is this?

`gbsa-grub` is an attempt to turn a familiar computational chemistry pattern

```text
ligand → docking → MD → GBSA
```

from a loose sequence of scripts into a reproducible workflow system.

The project combines three layers:

| layer | role |
|------|------|
| **gbsa-pipeline** | scientific modules (Docking, MD, GBSA) |
| **GRUBICY** | orchestration (actions, dependencies, execution order) |
| **signac** | job/data model |

Together, these layers make the workflow more than “something that can be run.” They make it something that can be **read**, **traced**, **reproduced**, and later **improved without collapsing under its own complexity**.

In that sense, `gbsa-grub` is not just about getting a docking pose or an MD trajectory. It is about building a workflow where every stage leaves behind a legible state and where downstream computation is driven by structure rather than manual intervention.

---

## why this exists

A lot of scientific pipelines begin the same way: one script solves one urgent problem, then a second script adapts the output, then a third script quietly assumes a filename convention, and eventually the whole workflow works mostly because the original author still remembers how it is supposed to work.

That is usually enough until the moment it is not.

What breaks first is often not the science, but the structure:
- files drift
- assumptions become implicit
- outputs stop being predictable
- reruns become hard to interpret
- extensions become risky

`gbsa-grub` exists to push back against that pattern. Its purpose is to make the workflow explicit enough that the repository itself explains the process.

---

## what is GRUBICY?

**GRUBICY** is the orchestration layer used here to define workflows as explicit actions with declared outputs, running on top of **signac**. In practice, it provides the structural glue: what can run, what depends on what, and which artifacts should exist after each step. It does not perform the chemistry itself; it organizes the execution of chemistry-facing modules.

That role matters more than it sounds. In projects like this, the hardest thing is often not writing a docking wrapper or launching a short MD job. The hard part is making those steps behave like parts of one system rather than unrelated local successes. GRUBICY is what gives the workflow a visible spine.

Within `gbsa-grub`, it is therefore the layer that turns “a set of scripts” into “a pipeline with declared behavior.”

---

## pipeline overview

```text
ligand
  ↓
docking
  ↓
md
  ↓
gbsa
```

The sequence itself is simple on purpose. The interesting part is not that these stages exist, but that they are now being encoded as inspectable transitions. Each stage should know:
- where its inputs come from
- what outputs it is expected to produce
- how those outputs will be consumed downstream
- whether it can be rerun without manual cleanup

That is the real design goal.

---

## current workflow status

| step | status | note |
|------|--------|------|
| **Docking** | wired | produces job-local docking artifacts |
| **MD** | wired | consumes docking output and writes trajectories |
| **GBSA** | wired | connected, but outputs are not yet good/reliable |

This means the first structural milestone has been reached: **data already moves forward through the workflow**. The repository is no longer just a plan. It now contains a functioning pipeline skeleton.

At the same time, the current version should still be read as an integration build. The fact that the chain runs is important, but it is not the same as saying the workflow is already scientifically mature.

---

## modules

### docking
- AutoDock Vina-based docking stage
- produces pose files and `result.json`
- currently wired and passing data forward
- docking box definition is still provisional in many cases

The docking stage already does what it needs to do for the current phase: it creates a standardized result that downstream steps can discover automatically. That may sound modest, but it is a decisive shift away from manual, one-off execution.

---

### md
- consumes protein structure + docking pose
- prefers **Meeko export** (`PDBQT → SDF`) over ad hoc ligand reconstruction
- parametrization via OpenMM + GAFF
- solvation and short MD
- writes trajectory artifacts such as `.xtc`

The MD step is where structural workflow design starts to matter much more. It has to bridge incompatible worlds: docking-like output on one side, simulation-ready input on the other. In many projects, that bridge is where reproducibility quietly dies. Here, the goal is to make that conversion path explicit and stable enough that it becomes part of the pipeline rather than hidden glue code.

---

### gbsa
- not yet fully wired
- planned target: `gmx_MMPBSA` / Uni-GBSA style endpoint scoring
- intended to consume MD outputs directly from the job structure

GBSA is already wired into the workflow, which means the structural connection exists and the stage can be executed as part of the pipeline. What is still missing is output quality: the current results are not yet in a state that should be treated as reliable endpoint scoring. Structurally, though, the important step has already happened — GBSA is no longer outside the pipeline.

---

## development philosophy

### phase 1 — wiring 
Everything must run end-to-end.  
Even imperfect or scientifically provisional outputs are acceptable as long as the structure is correct and the workflow behaves coherently.

### phase 2 — refinement (current)
Once the pipeline shape is stable:
- docking definitions can be improved
- parametrization details can be corrected
- force-field choices can be validated
- physical and numerical realism can be assessed properly

This separation is intentional. Structural bugs and scientific bugs become much harder to diagnose when they are mixed together. The current branch therefore optimizes first for **connectivity**, because a scientifically perfect module still has little value if it cannot live reliably inside a larger workflow.

---

## versioning

| version | meaning |
|--------|--------|
| 0.0.1 | pipeline wired |
| 0.1.0 | docking validated |
| 0.2.0 | MD validated |
| 0.3.0 | GBSA working |
| 0.4.0 | stabilization |
| 1.0.0 | scientifically usable |

This versioning scheme reflects the actual development logic of the project. A workflow should first become structurally real, then scientifically reliable.

---

## design principles

- **explicit > implicit**  
- **structure > scripts**  
- **reproducibility > convenience**  
- **modularity > monolith**

These are not decorative slogans. They are implementation rules. Whenever a shortcut hides state, makes dependencies harder to inspect, or forces downstream steps to guess what happened upstream, it is probably the wrong shortcut for this repository.

---

## repository layout

```text
gbsa-grub/
├── actions/
│   ├── 00_docking.py
│   ├── 01_md.py
│   └── 02_gbsa.py
├── pipeline.toml
├── workspace/
└── data/
```

The layout is centered around **actions** and **job-local artifacts**. The idea is that both code and outputs should reflect the same workflow structure.

---

## current limitations

- docking box selection is not yet robust
- MD is structurally wired but not scientifically validated
- GBSA is not yet fully connected
- the workflow is usable as an integration prototype, not yet as a production method

That is an acceptable state for this branch. The point of `wired` is to establish control over the workflow before claiming maturity.

---

## tl;dr

before:
```text
run.sh → ??? → results
```

after:
```text
jobs → actions → dependencies → outputs
```

`gbsa-grub` is the point where the workflow stops being a collection of remembered steps and starts becoming a system.
