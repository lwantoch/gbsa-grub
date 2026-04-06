# gbsa-grub

<p align="center">
  <b>structured molecular workflows</b><br>
  <sub>Docking → MD → GBSA</sub>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/status-dev-orange">
  <img src="https://img.shields.io/badge/version-0.1.0--dev-blue">
  <img src="https://img.shields.io/badge/core-gbsa--pipeline-black">
  <img src="https://img.shields.io/badge/orchestration-GRUBICY-purple">
</p>

---

## what is this?

`gbsa-grub` combines:

- **gbsa-pipeline** → scientific computation  
- **GRUBICY** → workflow orchestration  
- **signac** → job storage  

The goal is not to replace existing tools, but to connect them in a way that makes the full workflow understandable and reproducible.  
Instead of running isolated scripts, each computation becomes part of a structured pipeline with clearly defined inputs, outputs, and dependencies.

---

## why GRUBICY matters here

GRUBICY is the core architectural layer of this project.

Without it, workflows typically look like this:

```
scripts → files → confusion
```

With GRUBICY:

```
jobs → dependencies → pipeline
```

This difference is not cosmetic — it fundamentally changes how results are generated, tracked, and reused.

GRUBICY extends signac by introducing:

- explicit **parent–child relationships**
- support for **multi-stage pipelines**
- a **single source of truth** for workflow definition
- safe **schema evolution via migrations**
- **typed parameters** using Pydantic

This allows us to move from “running things” to actually **building a system**.

---

## core idea (GRUBICY)

```
Docking job
    ↓
MD job
    ↓
GBSA job
```

Each job:

- knows where it came from  
- stores its own results  
- can be reproduced independently  

This solves a key issue in scientific workflows:

```
same parameters + different parent ≠ same result
```

GRUBICY ensures these cases are handled correctly.

---

## pipeline definition

Instead of scattering logic across scripts, the pipeline is defined centrally:

```
pipeline.toml
```

This file describes:

- which stages exist  
- how they depend on each other  
- which parameters are required  

As a result, the workflow becomes inspectable and easier to evolve.

---

## current implementation

```
Docking  ✔ implemented
MD       ✘ planned
GBSA     ✘ planned
```

Current focus is on:

- stabilizing docking  
- defining clean interfaces for downstream stages  
- ensuring reproducibility from the beginning  

---

## docking step (implemented)

```
ligand → RDKit → Meeko → PDBQT
receptor → OpenBabel → PDBQT
→ Vina → pose
```

In more detail:

- ligands are normalized via RDKit  
- Meeko generates Vina-compatible PDBQT  
- receptors are converted using OpenBabel  
- Vina performs docking  
- outputs are stored in a structured directory  

All transformations are explicit and logged.

---

## job layout

```
job/
└── docking/
    ├── ligand_input.pdbqt
    ├── receptor.pdbqt
    ├── ligand_vina_out.pdbqt
    ├── vina.log
    └── result.json
```

Each job is self-contained, which makes debugging and reuse straightforward.

---

## result.json (minimal by design)

```
{
  "status": "success",
  "final_pose_pdbqt": "...",
  "best_score": -3.8,
  "log_file": "vina.log"
}
```

The JSON intentionally contains only essential information.

Everything else is stored in logs.

---

## logs

All detailed output is written to:

```
vina.log
```

This includes:

- full stdout  
- full stderr  
- scoring tables  
- warnings  

This separation keeps summaries clean while preserving full traceability.

---

## design principles

- minimal JSON output  
- full log preservation  
- explicit transformations  
- no hidden state  
- reproducible job structure  

---

## stack

- RDKit  
- Meeko  
- Open Babel  
- Vina  
- signac  
- GRUBICY  
- gbsa-pipeline  

---

## version

```
0.1.0-dev
```

This indicates:

- active development  
- unstable interfaces  
- focus on iteration  

---

## roadmap

```
0.2 → MD integration
0.3 → GBSA integration
1.0 → stable pipeline
```

---

## goal

```
ligand → docking → MD → trajectory → GBSA
```

The final system should allow:

- full traceability  
- easy reruns  
- structured analysis  

---

## tl;dr

before:

```
run.sh → ??? → results
```

after:

```
jobs → dependencies → reproducible pipeline
```

The shift is simple:

from running commands  
to building workflows
