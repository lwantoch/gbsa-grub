"""Parametrize and solvate the docked protein-ligand system in one action.

This action combines parametrization and solvation on purpose. The current
``gbsa_pipeline.solvation_openmm.solvate_openmm`` helper needs the in-memory
``forcefield`` and ``parmed_structure`` objects returned by ``parametrize``.
Those Python objects do not survive separate Grubicy actions, so splitting
parametrization and solvation would require a stable reload/cache interface
first. This runner therefore performs both steps in one process and writes the
stable ``setup/solvated.gro`` and ``setup/solvated.top`` files that downstream
MD actions can load independently.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from gbsa_pipeline.parametrization import ParametrizationInput, parametrize
from gbsa_pipeline.solvation_box import SolvationParams
from gbsa_pipeline.solvation_openmm import solvate_openmm
from grubicy import open_job_from_directory, parent_file
from grubicy.typed import load_action_params
from pydantic import BaseModel, ConfigDict

PARENT_MANIFEST = "export/result.json"


class SetupSystemParams(BaseModel):
    """Validated parameters for the setup action.

    Grubicy stores action parameters in the signac job statepoint, and
    ``load_action_params`` validates those values against this model at runtime.
    The setup action still reads molecular input files from the parent
    ``export/result.json`` manifest because those files are products of the
    previous action, not setup parameters. This model covers only the solvation
    options that belong to the current action. Extra parameters are forbidden so
    accidental statepoint/schema drift fails at the workflow boundary.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    water_model: str
    box_shape: str
    padding_nm: float
    ion_concentration: float | None = None
    neutralize: bool


def _read_json(path: str | Path) -> dict[str, Any]:
    """Read a JSON object from a stage manifest.

    The setup action consumes the export-stage manifest instead of reconstructing
    upstream paths manually. The decoded JSON value must be an object because the
    runner expects named fields such as ``protein_pdb`` and ``docked_sdf``. This
    helper performs only generic JSON validation while action-specific required
    keys are checked separately. Molecular systems are kept in normal files and
    are not embedded in the manifest.
    """
    manifest_path = Path(path)

    if not manifest_path.exists():
        raise FileNotFoundError(f"JSON manifest not found: {manifest_path}")

    if not manifest_path.is_file():
        raise ValueError(f"JSON manifest path is not a file: {manifest_path}")

    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in manifest: {manifest_path}") from exc

    if not isinstance(data, dict):
        raise ValueError(f"Expected JSON object in manifest: {manifest_path}")

    return data


def _write_json(path: Path, data: dict[str, Any]) -> Path:
    """Write a small JSON manifest with deterministic formatting.

    The setup action writes ``setup/result.json`` as the contract consumed by
    the first MD action. The manifest records the stable solvated GROMACS files,
    the intermediate parametrization files, and the input files used for this
    setup step. Serialization errors are not hidden because non-JSON-compatible
    manifest data is a workflow bug. Large molecular artefacts are written as
    separate files, not embedded into this JSON file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _read_parent_manifest(job: Any) -> dict[str, Any]:
    """Read ``export/result.json`` from the direct Grubicy parent job.

    Grubicy owns parent-job resolution through the dependency declared in
    ``pipeline.toml``. This helper asks Grubicy for the parent manifest file
    instead of reconstructing parent paths manually from signac details. JSON
    parsing remains local because ``export/result.json`` is a gbsa-grub stage
    manifest, not a Grubicy action-parameter model. Missing parent products fail
    before parametrization starts.
    """
    manifest_path = parent_file(job, PARENT_MANIFEST, must_exist=True)
    return _read_json(manifest_path)


def _required_path(data: dict[str, Any], key: str) -> Path:
    """Return a required manifest path as a resolved file path.

    The setup action needs the original protein PDB and the exported docked SDF.
    This helper keeps required-key, value-type, file-existence, and file-kind
    checks explicit at the workflow boundary. Missing keys and missing files
    fail before parametrization or solvation starts, giving a direct path to
    inspect. The function does not parse or validate molecular contents.
    """
    if key not in data:
        available = ", ".join(sorted(data)) or "<none>"
        raise KeyError(
            f"Parent manifest is missing required key: {key}\n"
            f"Required parent manifest: {PARENT_MANIFEST}\n"
            f"Available keys: {available}"
        )

    value = data[key]
    if not isinstance(value, str) or not value.strip():
        raise ValueError(
            f"Parent manifest key '{key}' must be a non-empty string path. "
            f"Got: {value!r}"
        )

    path = Path(value).resolve()

    if not path.exists():
        raise FileNotFoundError(f"Required parent file not found for '{key}': {path}")

    if not path.is_file():
        raise ValueError(f"Required parent path for '{key}' is not a file: {path}")

    return path


def main(directory: str | Path) -> None:
    """Parametrize, solvate, and write stable GROMACS setup products.

    The action opens the current job from the workspace directory supplied by
    row/Grubicy. Parent molecular inputs are read from the direct parent
    ``export/result.json`` manifest, while current solvation options are loaded
    through Grubicy's typed runtime helper. Parametrization and solvation remain
    delegated to ``gbsa_pipeline`` so this runner only coordinates workflow
    files, output directories, and the setup-stage manifest. The resulting
    ``gro_file`` and ``top_file`` keys are written for the first MD action.
    """
    job = open_job_from_directory(str(directory))
    params = load_action_params(job, SetupSystemParams)
    parent = _read_parent_manifest(job)

    protein_pdb = _required_path(parent, "protein_pdb")
    docked_sdf = _required_path(parent, "docked_sdf")

    setup_dir = Path(job.fn("setup"))
    parametrization_dir = setup_dir / "parametrization"
    solvation_dir = setup_dir / "solvation"

    setup_dir.mkdir(parents=True, exist_ok=True)
    parametrization_dir.mkdir(parents=True, exist_ok=True)
    solvation_dir.mkdir(parents=True, exist_ok=True)

    parametrized = parametrize(
        ParametrizationInput(
            protein_pdb=protein_pdb,
            ligand_sdf=docked_sdf,
            work_dir=parametrization_dir,
        )
    )

    solvated = solvate_openmm(
        parametrized=parametrized,
        params=SolvationParams(
            water_model=params.water_model,
            shape=params.box_shape,
            padding=params.padding_nm,
            ion_concentration=params.ion_concentration,
            neutralize=params.neutralize,
        ),
        output_gro=setup_dir / "solvated.gro",
        output_top=setup_dir / "solvated.top",
    )

    solvated_gro = Path(solvated.gro_file).resolve()
    solvated_top = Path(solvated.top_file).resolve()
    parametrized_gro = Path(parametrized.gro_file).resolve()
    parametrized_top = Path(parametrized.top_file).resolve()

    crystal_waters_pdb = (
        str(Path(parametrized.crystal_waters_pdb).resolve())
        if parametrized.crystal_waters_pdb is not None
        else None
    )

    _write_json(
        setup_dir / "result.json",
        {
            "protein_pdb": str(protein_pdb),
            "docked_sdf": str(docked_sdf),
            "solvated_gro": str(solvated_gro),
            "solvated_top": str(solvated_top),
            "gro_file": str(solvated_gro),
            "top_file": str(solvated_top),
            "parametrized_gro": str(parametrized_gro),
            "parametrized_top": str(parametrized_top),
            "crystal_waters_pdb": crystal_waters_pdb,
            "setup_dir": str(setup_dir),
            "parametrization_dir": str(parametrization_dir),
            "solvation_dir": str(solvation_dir),
            "water_model": params.water_model,
            "box_shape": params.box_shape,
            "padding_nm": params.padding_nm,
            "ion_concentration": params.ion_concentration,
            "neutralize": params.neutralize,
            "score": parent.get("score"),
            "rank": parent.get("rank"),
            "engine": parent.get("engine"),
        },
    )


if __name__ == "__main__":
    main(sys.argv[1])
