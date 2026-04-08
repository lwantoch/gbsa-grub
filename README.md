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

```
STATUS: pipeline wiring in progress
PRIORITY: connectivity > correctness
```

- full pipeline path is being connected end-to-end  
- modules are integrated, not yet validated  
- structure and reproducibility are the focus  

---

## what is this?

`gbsa-grub` connects three layers:

| layer | role |
|------|------|
| gbsa-pipeline | scientific logic |
| GRUBICY | workflow orchestration |
| signac | job/data model |

→ result: structured, reproducible molecular workflows

---

## why GRUBICY

without:

```
run.sh → scripts → files → chaos
```

with:

```
jobs → dependencies → pipeline
```

GRUBICY provides:

- explicit parent-child relationships  
- deterministic pipelines  
- central workflow definition (`pipeline.toml`)  
- reproducible execution  

---

## core idea

```
Docking → MD → GBSA
```

Each step:

- knows its parent  
- produces its own outputs  
- can be rerun independently  

---

## current status (wired)

| step     | status |
|----------|--------|
| Docking  | ✔ (wired, box arbitrary) |
| MD       | ✔ (wired) |
| GBSA     | ✘ |

---

## pipeline view

```
ligand
  ↓
docking
  ↓
md
  ↓
gbsa
```

---

## development strategy

### phase 1 — wiring (current)

- connect all modules  
- ensure data flows through pipeline  
- ignore scientific correctness  

### phase 2 — refinement

- fix chemistry  
- fix parametrization  
- validate results  

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

```
ligand → docking → MD → GBSA
```

---

## tl;dr

before:
```
run.sh → ??? → results
```

after:
```
jobs → dependencies → pipeline
```
