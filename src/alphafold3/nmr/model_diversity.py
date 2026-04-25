# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0.

"""AF3 ranked-model diversity metrics from explicit atom identities."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
import dataclasses
import math
from typing import Any

import numpy as np

from alphafold3.nmr import geometry_bundle


DIVERSITY_NEAR_DUPLICATE = "near_duplicate_models"
DIVERSITY_MODEST_LOCAL = "modest_local_diversity"
DIVERSITY_ENSEMBLE_LIKE = "ensemble_like"
DIVERSITY_NOT_APPLICABLE = "not_applicable"

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})
AROMATIC_RING_ATOMS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HID": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIE": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIP": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
METHYL_ATOMS = {
    "ALA": ("CB",),
    "VAL": ("CG1", "CG2"),
    "LEU": ("CD1", "CD2"),
    "ILE": ("CD1", "CG2"),
    "THR": ("CG2",),
    "MET": ("CE",),
}
AMIDE_HYDROGEN_NAMES = ("H", "HN")

MIN_DIVERSITY_MODELS = 2
MIN_SHARED_BACKBONE_ATOMS = 20
MIN_SHARED_RESIDUES = 5


@dataclasses.dataclass(frozen=True, slots=True)
class DiversityResult:
  pairwise_rows: list[dict[str, Any]]
  residue_variance_rows: list[dict[str, Any]]
  local_class_rows: list[dict[str, Any]]
  summary_rows: list[dict[str, Any]]
  diversity_grade: str
  median_pairwise_backbone_rmsd: float | None
  max_pairwise_backbone_rmsd: float | None
  p95_per_residue_variance: float | None
  shared_backbone_atoms: int
  shared_residues: int


def compute_model_diversity(
    atom_rows: Sequence[Mapping[str, Any]],
) -> DiversityResult:
  """Computes AF3-native diversity metrics across ranked models."""
  by_model = _atoms_by_model(atom_rows)
  model_ranks = sorted(by_model)
  if len(model_ranks) < MIN_DIVERSITY_MODELS:
    return _not_applicable_result(len(model_ranks), 0, 0)

  pairwise_rows = _pairwise_backbone_rmsd_rows(by_model)
  max_shared_backbone = max(
      (int(row["shared_backbone_atoms"]) for row in pairwise_rows), default=0
  )
  max_shared_residues = max(
      (int(row["shared_residues"]) for row in pairwise_rows), default=0
  )
  residue_variance_rows = _per_residue_variance_rows(by_model)
  if (
      max_shared_backbone < MIN_SHARED_BACKBONE_ATOMS
      or max_shared_residues < MIN_SHARED_RESIDUES
  ):
    return DiversityResult(
        pairwise_rows=pairwise_rows,
        residue_variance_rows=residue_variance_rows,
        local_class_rows=_local_class_rows(by_model, residue_variance_rows),
        summary_rows=[
            _summary_row(
                atom_rows=atom_rows,
                diversity_grade=DIVERSITY_NOT_APPLICABLE,
                mean_rmsd=None,
                median_rmsd=None,
                max_rmsd=None,
                median_variance=None,
                p95_variance=None,
                shared_atom_fraction=None,
                shared_backbone_atoms=max_shared_backbone,
                shared_residues=max_shared_residues,
            )
        ],
        diversity_grade=DIVERSITY_NOT_APPLICABLE,
        median_pairwise_backbone_rmsd=None,
        max_pairwise_backbone_rmsd=None,
        p95_per_residue_variance=None,
        shared_backbone_atoms=max_shared_backbone,
        shared_residues=max_shared_residues,
    )

  rmsds = [
      float(row["backbone_rmsd"])
      for row in pairwise_rows
      if row["backbone_rmsd"] != ""
  ]
  shared_fractions = [
      float(row["shared_atom_fraction"])
      for row in pairwise_rows
      if row.get("shared_atom_fraction") not in ("", None)
  ]
  variances = [
      float(row["coordinate_variance"])
      for row in residue_variance_rows
      if row["coordinate_variance"] != ""
  ]
  mean_rmsd = float(np.mean(rmsds)) if rmsds else None
  median_rmsd = float(np.median(rmsds)) if rmsds else None
  max_rmsd = float(np.max(rmsds)) if rmsds else None
  median_variance = float(np.median(variances)) if variances else None
  p95_variance = float(np.percentile(variances, 95)) if variances else None
  mean_shared_fraction = float(np.mean(shared_fractions)) if shared_fractions else None
  grade = diversity_grade(median_rmsd, p95_variance)
  local_class_rows = _local_class_rows(by_model, residue_variance_rows)
  return DiversityResult(
      pairwise_rows=pairwise_rows,
      residue_variance_rows=residue_variance_rows,
      local_class_rows=local_class_rows,
      summary_rows=[
          _summary_row(
              atom_rows=atom_rows,
              diversity_grade=grade,
              mean_rmsd=mean_rmsd,
              median_rmsd=median_rmsd,
              max_rmsd=max_rmsd,
              median_variance=median_variance,
              p95_variance=p95_variance,
              shared_atom_fraction=mean_shared_fraction,
              shared_backbone_atoms=max_shared_backbone,
              shared_residues=max_shared_residues,
          )
      ],
      diversity_grade=grade,
      median_pairwise_backbone_rmsd=median_rmsd,
      max_pairwise_backbone_rmsd=max_rmsd,
      p95_per_residue_variance=p95_variance,
      shared_backbone_atoms=max_shared_backbone,
      shared_residues=max_shared_residues,
  )


def diversity_grade(
    median_pairwise_backbone_rmsd: float | None,
    p95_per_residue_variance: float | None,
) -> str:
  if median_pairwise_backbone_rmsd is None or p95_per_residue_variance is None:
    return DIVERSITY_NOT_APPLICABLE
  if (
      median_pairwise_backbone_rmsd >= 0.5
      or p95_per_residue_variance >= 0.25
  ):
    return DIVERSITY_ENSEMBLE_LIKE
  if (
      median_pairwise_backbone_rmsd < 0.1
      and p95_per_residue_variance < 0.05
  ):
    return DIVERSITY_NEAR_DUPLICATE
  return DIVERSITY_MODEST_LOCAL


def _not_applicable_result(
    n_models: int, shared_backbone_atoms: int, shared_residues: int
) -> DiversityResult:
  row = {
      "system_id": "",
      "n_models": n_models,
      "mean_pairwise_backbone_rmsd": "",
      "median_pairwise_backbone_rmsd": "",
      "max_pairwise_backbone_rmsd": "",
      "median_per_residue_variance": "",
      "p95_per_residue_variance": "",
      "shared_atom_fraction": "",
      "shared_backbone_atoms": shared_backbone_atoms,
      "shared_residues": shared_residues,
      "diversity_grade": DIVERSITY_NOT_APPLICABLE,
  }
  return DiversityResult(
      pairwise_rows=[],
      residue_variance_rows=[],
      local_class_rows=[],
      summary_rows=[row],
      diversity_grade=DIVERSITY_NOT_APPLICABLE,
      median_pairwise_backbone_rmsd=None,
      max_pairwise_backbone_rmsd=None,
      p95_per_residue_variance=None,
      shared_backbone_atoms=shared_backbone_atoms,
      shared_residues=shared_residues,
  )


def _atoms_by_model(
    atom_rows: Sequence[Mapping[str, Any]],
) -> dict[str, dict[tuple[str, ...], Mapping[str, Any]]]:
  by_model: dict[str, dict[tuple[str, ...], Mapping[str, Any]]] = defaultdict(dict)
  for row in atom_rows:
    by_model[str(row["model_rank"])][geometry_bundle.identity_without_model(row)] = row
  return dict(by_model)


def _pairwise_backbone_rmsd_rows(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> list[dict[str, Any]]:
  model_ranks = sorted(by_model)
  pairs = _model_pairs(model_ranks)
  rows: list[dict[str, Any]] = []
  for model_a, model_b in pairs:
    atoms_a = by_model[model_a]
    atoms_b = by_model[model_b]
    shared_keys = sorted(
        key for key in (set(atoms_a) & set(atoms_b)) if key[4] in BACKBONE_ATOMS
    )
    shared_residues = {
        (key[1], key[2], key[3]) for key in shared_keys
    }
    shared_atom_fraction = (
        len(shared_keys) / min(len(atoms_a), len(atoms_b))
        if min(len(atoms_a), len(atoms_b))
        else 0.0
    )
    if len(shared_keys) >= 3:
      coords_a = np.array([_coords(atoms_a[key]) for key in shared_keys])
      coords_b = np.array([_coords(atoms_b[key]) for key in shared_keys])
      rmsd = _aligned_rmsd(coords_a, coords_b)
    else:
      rmsd = None
    rows.append(
        {
            "system_id": _first_system_id(atoms_a.values()),
            "model_rank_a": model_a,
            "model_rank_b": model_b,
            "shared_backbone_atoms": len(shared_keys),
            "shared_residues": len(shared_residues),
            "shared_atom_fraction": shared_atom_fraction,
            "backbone_rmsd": rmsd if rmsd is not None else "",
        }
    )
  return rows


def _model_pairs(model_ranks: Sequence[str]) -> list[tuple[str, str]]:
  if len(model_ranks) <= 10:
    return [
        (model_a, model_b)
        for i, model_a in enumerate(model_ranks)
        for model_b in model_ranks[i + 1 :]
    ]
  rank_1_pairs = [(model_ranks[0], model) for model in model_ranks[1:]]
  adjacent_pairs = list(zip(model_ranks, model_ranks[1:]))
  return sorted(set(rank_1_pairs + adjacent_pairs))


def _per_residue_variance_rows(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> list[dict[str, Any]]:
  model_ranks = sorted(by_model)
  atom_keys = sorted(set().union(*(set(rows) for rows in by_model.values())))
  per_residue: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
  per_residue_backbone: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
  per_residue_sidechain: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
  for atom_key in atom_keys:
    present = [by_model[model].get(atom_key) for model in model_ranks]
    if any(row is None for row in present):
      continue
    rows = [row for row in present if row is not None]
    coords = np.array([_coords(row) for row in rows])
    variance = _coordinate_variance(coords)
    residue_key = (
        str(rows[0]["system_id"]),
        atom_key[1],
        atom_key[2],
        atom_key[3],
        str(rows[0].get("residue_name", "")),
    )
    per_residue[residue_key].append(variance)
    if atom_key[4] in BACKBONE_ATOMS:
      per_residue_backbone[residue_key].append(variance)
    else:
      per_residue_sidechain[residue_key].append(variance)

  rows: list[dict[str, Any]] = []
  for residue_key in sorted(per_residue):
    values = per_residue[residue_key]
    backbone_values = per_residue_backbone.get(residue_key, [])
    sidechain_values = per_residue_sidechain.get(residue_key, [])
    rows.append(
        {
            "system_id": residue_key[0],
            "chain_id": residue_key[1],
            "residue_number": residue_key[2],
            "insertion_code": residue_key[3],
            "residue_name": residue_key[4],
            "coordinate_variance": float(np.mean(values)),
            "backbone_variance": (
                float(np.mean(backbone_values)) if backbone_values else ""
            ),
            "sidechain_variance": (
                float(np.mean(sidechain_values)) if sidechain_values else ""
            ),
            "shared_atom_count": len(values),
            "n_models": len(model_ranks),
        }
    )
  return rows


def _local_class_rows(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
    residue_variance_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  rows = [
      _variance_class_summary("backbone", residue_variance_rows, "backbone_variance"),
      _variance_class_summary("sidechain", residue_variance_rows, "sidechain_variance"),
      _aromatic_summary(by_model),
      _methyl_summary(by_model),
      _amide_summary(by_model),
  ]
  return rows


def _variance_class_summary(
    class_name: str,
    residue_variance_rows: Sequence[Mapping[str, Any]],
    column: str,
) -> dict[str, Any]:
  values = [
      float(row[column])
      for row in residue_variance_rows
      if row.get(column) not in ("", None)
  ]
  return _class_summary_from_values(class_name, values)


def _aromatic_summary(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> dict[str, Any]:
  centroid_variances: list[float] = []
  plane_spreads: list[float] = []
  for residue_atoms in _residue_atom_groups(by_model):
    residue_name = residue_atoms.residue_name
    ring_atoms = AROMATIC_RING_ATOMS.get(residue_name)
    if not ring_atoms:
      continue
    coords_by_model = []
    normals = []
    for model_rank in residue_atoms.model_ranks:
      atoms = residue_atoms.atoms_by_model[model_rank]
      if not all(atom in atoms for atom in ring_atoms):
        break
      coords = np.array([_coords(atoms[atom]) for atom in ring_atoms])
      coords_by_model.append(coords.mean(axis=0))
      normals.append(_plane_normal(coords))
    else:
      centroid_variances.append(_coordinate_variance(np.array(coords_by_model)))
      plane_spreads.append(_angle_spread(normals))
  row = _class_summary_from_values("aromatic_ring_centroid", centroid_variances)
  row["plane_angle_mean_degrees"] = (
      float(np.mean(plane_spreads)) if plane_spreads else ""
  )
  row["plane_angle_p95_degrees"] = (
      float(np.percentile(plane_spreads, 95)) if plane_spreads else ""
  )
  return row


def _methyl_summary(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> dict[str, Any]:
  values: list[float] = []
  for residue_atoms in _residue_atom_groups(by_model):
    methyl_atoms = METHYL_ATOMS.get(residue_atoms.residue_name)
    if not methyl_atoms:
      continue
    for atom_name in methyl_atoms:
      coords = []
      for model_rank in residue_atoms.model_ranks:
        atom = residue_atoms.atoms_by_model[model_rank].get(atom_name)
        if atom is None:
          break
        coords.append(_coords(atom))
      else:
        values.append(_coordinate_variance(np.array(coords)))
  return _class_summary_from_values("methyl_carbon", values)


def _amide_summary(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> dict[str, Any]:
  angle_spreads: list[float] = []
  for residue_atoms in _residue_atom_groups(by_model):
    vectors = []
    for model_rank in residue_atoms.model_ranks:
      atoms = residue_atoms.atoms_by_model[model_rank]
      hydrogen = next((atoms[name] for name in AMIDE_HYDROGEN_NAMES if name in atoms), None)
      nitrogen = atoms.get("N")
      if hydrogen is None or nitrogen is None:
        break
      vector = _coords(hydrogen) - _coords(nitrogen)
      norm = np.linalg.norm(vector)
      if norm == 0:
        break
      vectors.append(vector / norm)
    else:
      angle_spreads.append(_angle_spread(vectors))
  row = _class_summary_from_values("amide_vector_angle", angle_spreads)
  row["units"] = "degrees"
  return row


def _class_summary_from_values(class_name: str, values: Sequence[float]) -> dict[str, Any]:
  if not values:
    return {
        "system_id": "",
        "class": class_name,
        "n_residues": 0,
        "mean_variance": "",
        "p95_variance": "",
        "max_variance": "",
        "high_diversity_fraction": "",
        "status": "not_applicable",
    }
  arr = np.array(values, dtype=float)
  return {
      "system_id": "",
      "class": class_name,
      "n_residues": int(arr.size),
      "mean_variance": float(np.mean(arr)),
      "p95_variance": float(np.percentile(arr, 95)),
      "max_variance": float(np.max(arr)),
      "high_diversity_fraction": float(np.mean(arr >= 0.25)),
      "status": "ok",
  }


@dataclasses.dataclass(frozen=True, slots=True)
class _ResidueAtomGroup:
  residue_name: str
  model_ranks: tuple[str, ...]
  atoms_by_model: Mapping[str, Mapping[str, Mapping[str, Any]]]


def _residue_atom_groups(
    by_model: Mapping[str, Mapping[tuple[str, ...], Mapping[str, Any]]],
) -> list[_ResidueAtomGroup]:
  model_ranks = tuple(sorted(by_model))
  grouped: dict[tuple[str, str, str, str], dict[str, dict[str, Mapping[str, Any]]]] = defaultdict(
      lambda: defaultdict(dict)
  )
  residue_names: dict[tuple[str, str, str, str], str] = {}
  for model_rank, atoms in by_model.items():
    for atom_key, row in atoms.items():
      residue_key = (atom_key[0], atom_key[1], atom_key[2], atom_key[3])
      grouped[residue_key][model_rank][atom_key[4]] = row
      residue_names[residue_key] = str(row.get("residue_name", ""))
  groups = []
  for residue_key in sorted(grouped):
    atoms_by_model = grouped[residue_key]
    if not all(model_rank in atoms_by_model for model_rank in model_ranks):
      continue
    groups.append(
        _ResidueAtomGroup(
            residue_name=residue_names.get(residue_key, ""),
            model_ranks=model_ranks,
            atoms_by_model=atoms_by_model,
        )
    )
  return groups


def _summary_row(
    *,
    atom_rows: Sequence[Mapping[str, Any]],
    diversity_grade: str,
    mean_rmsd: float | None,
    median_rmsd: float | None,
    max_rmsd: float | None,
    median_variance: float | None,
    p95_variance: float | None,
    shared_atom_fraction: float | None,
    shared_backbone_atoms: int,
    shared_residues: int,
) -> dict[str, Any]:
  model_ranks = sorted({str(row["model_rank"]) for row in atom_rows})
  system_id = str(atom_rows[0]["system_id"]) if atom_rows else ""
  return {
      "system_id": system_id,
      "n_models": len(model_ranks),
      "mean_pairwise_backbone_rmsd": (
          mean_rmsd if mean_rmsd is not None else ""
      ),
      "median_pairwise_backbone_rmsd": median_rmsd if median_rmsd is not None else "",
      "max_pairwise_backbone_rmsd": max_rmsd if max_rmsd is not None else "",
      "median_per_residue_variance": (
          median_variance if median_variance is not None else ""
      ),
      "p95_per_residue_variance": (
          p95_variance if p95_variance is not None else ""
      ),
      "shared_atom_fraction": shared_atom_fraction,
      "shared_backbone_atoms": shared_backbone_atoms,
      "shared_residues": shared_residues,
      "diversity_grade": diversity_grade,
  }


def _coords(row: Mapping[str, Any]) -> np.ndarray:
  return np.array([float(row["x"]), float(row["y"]), float(row["z"])], dtype=float)


def _coordinate_variance(coords: np.ndarray) -> float:
  centroid = coords.mean(axis=0)
  return float(np.mean(np.sum((coords - centroid) ** 2, axis=1)))


def _aligned_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
  if coords_a.shape != coords_b.shape:
    raise ValueError("coordinate arrays must have matching shape")
  center_a = coords_a.mean(axis=0)
  center_b = coords_b.mean(axis=0)
  a = coords_a - center_a
  b = coords_b - center_b
  covariance = a.T @ b / coords_a.shape[0]
  u, _, vt = np.linalg.svd(covariance)
  rotation = vt.T @ u.T
  if np.linalg.det(rotation) < 0:
    vt[-1, :] *= -1
    rotation = vt.T @ u.T
  aligned = a @ rotation + center_b
  return float(np.sqrt(np.mean(np.sum((aligned - coords_b) ** 2, axis=1))))


def _plane_normal(coords: np.ndarray) -> np.ndarray:
  centered = coords - coords.mean(axis=0)
  _, _, vh = np.linalg.svd(centered)
  normal = vh[-1]
  norm = np.linalg.norm(normal)
  return normal / norm if norm else normal


def _angle_spread(vectors: Sequence[np.ndarray]) -> float:
  if len(vectors) < 2:
    return 0.0
  base = vectors[0]
  angles = []
  for vector in vectors[1:]:
    dot = abs(float(np.dot(base, vector)))
    dot = min(1.0, max(-1.0, dot))
    angles.append(math.degrees(math.acos(dot)))
  return float(np.max(angles)) if angles else 0.0


def _first_system_id(rows: Sequence[Mapping[str, Any]] | Any) -> str:
  for row in rows:
    return str(row["system_id"])
  return ""
