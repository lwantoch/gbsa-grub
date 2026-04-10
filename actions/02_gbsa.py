# /home/grheco/PycharmProjects/gbsa-grub/actions/02_gbsa.py

"""GRUBICY action for wiring a GBSA step on top of an MD parent job.

Design goals (wired phase):
- connectivity > correctness
- tolerant path resolution from MD outputs
- automatic group detection from index.ndx (with override)
- deterministic artifacts:
    - gbsa/mmpbsa.in
    - gbsa/gmx_mmpbsa.stdout.log
    - gbsa/gmx_mmpbsa.stderr.log
    - gbsa/result.json
    - gbsa_action/result.json
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

from gbsa_pipeline.unigbsa import (
    GBParams,
    GeneralParams,
    GmxMmpbsaRequest,
    MMPBSAConfig,
    PBParams,
    run_gmx_mmpbsa,
)

STATEPOINT_FILENAME = "signac_statepoint.json"


# -----------------------------------------------------------------------------
# IO helpers
# -----------------------------------------------------------------------------


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required {label} file does not exist: {path}")
    return path


def _resolve_relative_to(base_dir: Path, value: str | None) -> Path | None:
    if not value:
        return None
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (base_dir / candidate).resolve()


def _first_existing(paths: list[Path | None], label: str) -> Path:
    for path in paths:
        if path is not None and path.exists():
            return path
    attempted = [str(p) for p in paths if p is not None]
    raise FileNotFoundError(f"Could not resolve {label}. Tried: {attempted}")


def _get_nested(mapping: dict[str, Any], *keys: str) -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return None
        current = current[key]
    return current


# -----------------------------------------------------------------------------
# NDX parsing / group detection
# -----------------------------------------------------------------------------


def _parse_ndx_groups(index_file: Path) -> list[str]:
    groups: list[str] = []
    for line in index_file.read_text().splitlines():
        line = line.strip()
        if line.startswith("[") and line.endswith("]"):
            groups.append(line[1:-1].strip())
    return groups


def _auto_detect_groups(index_file: Path) -> tuple[int, int]:
    """Heuristic detection of receptor and ligand groups from index.ndx."""
    groups = _parse_ndx_groups(index_file)

    receptor_idx = None
    ligand_idx = None

    for i, name in enumerate(groups):
        lname = name.lower()

        if receptor_idx is None and "protein" in lname:
            receptor_idx = i

        if ligand_idx is None and any(x in lname for x in ["lig", "ligand", "mol"]):
            ligand_idx = i

    # fallbacks
    if receptor_idx is None:
        receptor_idx = 0

    if ligand_idx is None:
        if len(groups) > 1:
            ligand_idx = 1
        else:
            raise RuntimeError("Could not detect ligand group in index.ndx")

    return receptor_idx, ligand_idx


# -----------------------------------------------------------------------------
# State / parent loading
# -----------------------------------------------------------------------------


def _load_statepoint(job_dir: Path) -> dict[str, Any]:
    return _read_json(_require_file(job_dir / STATEPOINT_FILENAME, "statepoint"))


def _load_parent_md_manifest(
    job_dir: Path, statepoint: dict[str, Any]
) -> tuple[Path, dict[str, Any]]:
    rel = statepoint.get("parent_md_result", "md/result.json")
    path = _resolve_relative_to(job_dir, rel)
    if path is None:
        raise ValueError("parent_md_result could not be resolved")
    path = _require_file(path, "parent MD result manifest")
    return path, _read_json(path)


# -----------------------------------------------------------------------------
# Input resolution
# -----------------------------------------------------------------------------


def _resolve_md_inputs(
    job_dir: Path,
    statepoint: dict[str, Any],
    parent_manifest: dict[str, Any],
) -> dict[str, Path]:
    """Resolve complex_structure, trajectory, topology, index_file."""

    explicit_complex = _resolve_relative_to(
        job_dir, statepoint.get("gbsa_complex_structure")
    )
    explicit_traj = _resolve_relative_to(job_dir, statepoint.get("gbsa_trajectory"))
    explicit_top = _resolve_relative_to(job_dir, statepoint.get("gbsa_topology"))
    explicit_ndx = _resolve_relative_to(job_dir, statepoint.get("gbsa_index_file"))

    md_dir = _resolve_relative_to(
        job_dir, _get_nested(parent_manifest, "md", "directory")
    )

    complex_structure = _first_existing(
        [
            explicit_complex,
            _resolve_relative_to(
                job_dir, _get_nested(parent_manifest, "md", "final_gro")
            ),
            _resolve_relative_to(
                job_dir, _get_nested(parent_manifest, "md", "final_pdb")
            ),
            md_dir / "md_final.gro" if md_dir else None,
            md_dir / "md_final.pdb" if md_dir else None,
            job_dir / "md" / "native_outputs" / "md.gro",
            job_dir / "md" / "native_outputs" / "md_final.gro",
        ],
        "complex structure",
    )

    trajectory = _first_existing(
        [
            explicit_traj,
            _resolve_relative_to(
                job_dir, _get_nested(parent_manifest, "md", "trajectory")
            ),
            md_dir / "trajectory.xtc" if md_dir else None,
            job_dir / "md" / "native_outputs" / "trajectory.xtc",
        ],
        "trajectory",
    )

    topology = _first_existing(
        [
            explicit_top,
            _resolve_relative_to(
                job_dir, _get_nested(parent_manifest, "md", "topology")
            ),
            md_dir / "topol.top" if md_dir else None,
            job_dir / "md" / "native_outputs" / "topol.top",
        ],
        "topology",
    )

    index_file = _first_existing(
        [
            explicit_ndx,
            _resolve_relative_to(
                job_dir, _get_nested(parent_manifest, "md", "index_file")
            ),
            md_dir / "index.ndx" if md_dir else None,
            job_dir / "md" / "native_outputs" / "index.ndx",
        ],
        "index file",
    )

    return {
        "complex_structure": _require_file(complex_structure, "complex"),
        "trajectory": _require_file(trajectory, "trajectory"),
        "topology": _require_file(topology, "topology"),
        "index_file": _require_file(index_file, "index"),
    }


# -----------------------------------------------------------------------------
# Config builder
# -----------------------------------------------------------------------------


def _build_config(sp: dict[str, Any]) -> MMPBSAConfig:
    general = GeneralParams(
        startframe=sp.get("gbsa_startframe", 1),
        endframe=sp.get("gbsa_endframe", 9_999_999),
        interval=sp.get("gbsa_interval", 1),
        temperature=sp.get("gbsa_temperature", 298.15),
        verbose=sp.get("gbsa_verbose", 1),
        extra=sp.get("gbsa_general_extra", {}),
    )

    gb = (
        GBParams(
            igb=sp.get("gbsa_igb", 5),
            extra=sp.get("gbsa_gb_extra", {}),
        )
        if sp.get("gbsa_use_gb", True)
        else None
    )

    pb = (
        PBParams(
            ipb=sp.get("gbsa_ipb", 2),
            extra=sp.get("gbsa_pb_extra", {}),
        )
        if sp.get("gbsa_use_pb", True)
        else None
    )

    return MMPBSAConfig(general=general, gb=gb, pb=pb)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> int:
    job_dir = Path.cwd()
    action_dir = job_dir / "gbsa_action"
    action_dir.mkdir(exist_ok=True)

    print("[GBSA] start")

    try:
        sp = _load_statepoint(job_dir)
        parent_path, parent = _load_parent_md_manifest(job_dir, sp)
        inputs = _resolve_md_inputs(job_dir, sp, parent)

        # --- AUTO GROUP DETECTION ---
        receptor_group, ligand_group = _auto_detect_groups(inputs["index_file"])

        # optional override
        if "receptor_group" in sp:
            receptor_group = int(sp["receptor_group"])
        if "ligand_group" in sp:
            ligand_group = int(sp["ligand_group"])

        config = _build_config(sp)

        request = GmxMmpbsaRequest(
            complex_structure=inputs["complex_structure"],
            trajectory=inputs["trajectory"],
            topology=inputs["topology"],
            index_file=inputs["index_file"],
            receptor_group=receptor_group,
            ligand_group=ligand_group,
            output_dir=job_dir / "gbsa",
            config=config,
            gmx_mmpbsa=sp.get("gmx_mmpbsa_binary", "gmx_MMPBSA"),
            extra_args=tuple(sp.get("gmx_mmpbsa_extra_args", [])),
        )

        result = run_gmx_mmpbsa(request)

        _write_json(
            action_dir / "result.json",
            {
                "status": "success" if result.ok else "failed",
                "receptor_group": receptor_group,
                "ligand_group": ligand_group,
                "inputs": {k: str(v) for k, v in inputs.items()},
                "gbsa_result": result.result_json,
                "returncode": result.returncode,
            },
        )

        print("[GBSA] done")
        return 0 if result.ok else 1

    except Exception as e:
        _write_json(
            action_dir / "result.json",
            {"status": "error", "error": str(e)},
        )
        print("[GBSA] ERROR:", e, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
