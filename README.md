# gbsa-grub

[![Status](https://img.shields.io/badge/status-under%20development-orange)](https://github.com/lwantoch/gbsa-grub)
[![Version](https://img.shields.io/badge/version-0.1.0--dev-blue)](https://github.com/lwantoch/gbsa-grub)
[![Python](https://img.shields.io/badge/python-3.12+-blue)](https://www.python.org/)
[![Workflow](https://img.shields.io/badge/workflow-GRUBICY%20%2B%20signac-6f42c1)](https://github.com/lwantoch/gbsa-grub)

---

## Status

Version: 0.1.0  
Status: under development  
Maturity: experimental  

---

## Overview

gbsa-grub is a workflow-oriented project that aims to reconstruct the full gbsa-pipeline inside a GRUBICY / signac execution model.

Instead of running a monolithic pipeline, each step is exposed as a structured action:

- explicit inputs
- explicit outputs
- reproducible job directories
- inspectable intermediate states

---

## Goal

The long-term goal is:

> represent the entire gbsa-pipeline as a modular, restartable, inspectable workflow

Planned stages:

- ligand preparation  
- docking  
- parametrization  
- solvation  
- minimization  
- MD  
- GBSA  

---

## Current Implementation

Currently implemented:

- docking action

The docking step:

- loads ligand via RDKit  
- converts ligand → PDBQT via Meeko  
- converts receptor → PDBQT via Open Babel  
- runs AutoDock Vina  
- stores results in structured JSON  
- writes full logs to file  

---

## Design Principles

1. Minimal JSON output  
   → workflow-friendly, stable, machine-readable  

2. Full logs on disk  
   → no bloated JSON  

3. Explicit intermediate files  
   → easy debugging  

4. Reuse gbsa-pipeline  
   → no duplicated chemistry logic  

5. Workflow-first architecture  
   → not scripts, but composable steps  

---

## Example Output

{
  "status": "success",
  "engine": "vina",
  "prepared_ligand_pdbqt": "...",
  "final_pose_pdbqt": "...",
  "best_score": -7.4,
  "log_file": "..."
}

---

## Versioning Strategy

While under development:

- use 0.x.y versions

Examples:

- 0.1.0 → initial docking integration  
- 0.2.0 → GBSA integration  
- 0.3.0 → MD integration  

Only move to 1.0.0 when:

- APIs are stable  
- workflow is complete  
- outputs are consistent  
- documentation is solid  

---

## Roadmap

Short-term:

- stabilize docking  
- define output schema  
- improve logging consistency  

Mid-term:

- integrate GBSA  
- integrate MD  
- unify pipeline logic  

Long-term:

- full gbsa-pipeline as GRUBICY workflow  

---

## Summary

gbsa-grub is an experimental workflow layer around gbsa-pipeline.

It is:

- modular  
- reproducible  
- inspectable  

Current state:

→ early development (0.1.0)

---

## Metadata

Project: gbsa-grub  
Version: 0.1.0  
Status: under development  
