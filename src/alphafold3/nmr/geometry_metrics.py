# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0.

"""Science-facing AF3-native geometry readiness metrics."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
import math
from typing import Any

import numpy as np

from alphafold3.nmr import geometry_bundle
from alphafold3.nmr import model_diversity


CONFIDENCE_ALIGNED = "aligned"
CONFIDENCE_WEAK = "weak"
CONFIDENCE_DECOUPLED = "decoupled"
CONFIDENCE_NOT_APPLICABLE = "not_applicable"
MIN_CONFIDENCE_CORRELATION_RESIDUES = 30


def bundle_readiness_summary(
    *,
    system_id: str,
    tables: geometry_bundle.BundleTables,
    diversity: model_diversity.DiversityResult,
    confidence_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  confidence_status = (
      confidence_rows[0]["confidence_geometry_agreement"]
      if confidence_rows
      else CONFIDENCE_NOT_APPLICABLE
  )
  return [
      {
          "system_id": system_id,
          "n_models": len({row["model_rank"] for row in tables.atom_rows}),
          "n_chains": len(
              {
                  (row["model_rank"], row["chain_id"])
                  for row in tables.chain_rows
              }
          ),
          "n_residues": len(tables.residue_rows),
          "n_atoms": len(tables.atom_rows),
          "readiness_verdict": tables.readiness_verdict,
          "identity_grade": tables.identity_grade,
          "geometry_grade": tables.geometry_grade,
          "diversity_grade": diversity.diversity_grade,
          "confidence_geometry_agreement": confidence_status,
          "blocked_reasons": ";".join(tables.blocked_reasons),
          "confidence_linkage_status": _confidence_linkage_status(
              tables.extraction_metrics
          ),
          "atom_table_hash": tables.atom_table_hash,
      }
  ]


def identity_metrics_rows(
    extraction_metrics: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  rows = []
  for row in sorted(
      extraction_metrics,
      key=lambda metric: (metric.get("system_id", ""), metric.get("model_rank", "")),
  ):
    rows.append(
        {
            "system_id": row.get("system_id", ""),
            "model_rank": row.get("model_rank", ""),
            "n_chains": row.get("chain_count", 0),
            "n_residues": row.get("residue_count", 0),
            "n_atoms": row.get("atom_count", 0),
            "duplicate_canonical_atom_keys": row.get(
                "duplicate_atom_identity_count", 0
            ),
            "auth_vs_label_chain_discrepancies": row.get(
                "auth_vs_label_chain_discrepancies", 0
            ),
            "auth_vs_label_residue_discrepancies": row.get(
                "auth_vs_label_residue_discrepancies", 0
            ),
            "required_atom_site_columns_present": row.get(
                "required_atom_site_columns_present", False
            ),
            "confidence_atom_count_matches": row.get(
                "confidence_atom_count_matches", False
            ),
        }
    )
  return rows


def geometry_completeness_rows(
    extraction_metrics: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  rows = []
  for row in sorted(
      extraction_metrics,
      key=lambda metric: (metric.get("system_id", ""), metric.get("model_rank", "")),
  ):
    rows.append(
        {
            "system_id": row.get("system_id", ""),
            "model_rank": row.get("model_rank", ""),
            "atom_count": row.get("atom_count", 0),
            "chain_count": row.get("chain_count", 0),
            "residue_count": row.get("residue_count", 0),
            "invalid_coordinate_count": row.get("invalid_coordinate_count", 0),
            "disallowed_altloc_count": row.get("disallowed_altloc_count", 0),
            "missing_atom_site_columns": row.get("missing_atom_site_columns", ""),
        }
    )
  return rows


def confidence_geometry_rows(
    atom_rows: Sequence[Mapping[str, Any]],
    residue_variance_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  plddt_by_residue = _mean_plddt_by_residue(atom_rows)
  paired = []
  for row in residue_variance_rows:
    key = (
        str(row["system_id"]),
        str(row["chain_id"]),
        str(row["residue_number"]),
        str(row["insertion_code"]),
        str(row["residue_name"]),
    )
    plddt = plddt_by_residue.get(key)
    variance = row.get("coordinate_variance")
    if plddt is None or variance in ("", None):
      continue
    paired.append((float(plddt), float(variance)))

  if len(paired) < MIN_CONFIDENCE_CORRELATION_RESIDUES:
    return [
        {
            "system_id": _first_system_id(atom_rows),
            "n_residues": len(paired),
            "spearman_plddt_vs_variance": "",
            "confidence_geometry_agreement": CONFIDENCE_NOT_APPLICABLE,
        }
    ]
  plddt_values = np.array([item[0] for item in paired], dtype=float)
  variance_values = np.array([item[1] for item in paired], dtype=float)
  rho = _spearman(plddt_values, variance_values)
  return [
      {
          "system_id": _first_system_id(atom_rows),
          "n_residues": len(paired),
          "spearman_plddt_vs_variance": rho,
          "confidence_geometry_agreement": confidence_geometry_agreement(rho),
      }
  ]


def confidence_geometry_agreement(rho: float | None) -> str:
  if rho is None or not math.isfinite(rho):
    return CONFIDENCE_NOT_APPLICABLE
  if rho <= -0.5:
    return CONFIDENCE_ALIGNED
  if abs(rho) < 0.2:
    return CONFIDENCE_DECOUPLED
  return CONFIDENCE_WEAK


def metrics_json(
    *,
    system_id: str,
    tables: geometry_bundle.BundleTables,
    diversity: model_diversity.DiversityResult,
    confidence_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
  confidence_status = (
      confidence_rows[0]["confidence_geometry_agreement"]
      if confidence_rows
      else CONFIDENCE_NOT_APPLICABLE
  )
  return {
      "atom_table_hash": tables.atom_table_hash,
      "blocked_reasons": list(tables.blocked_reasons),
      "confidence_geometry_agreement": confidence_status,
      "diversity_grade": diversity.diversity_grade,
      "geometry_grade": tables.geometry_grade,
      "identity_grade": tables.identity_grade,
      "readiness_verdict": tables.readiness_verdict,
      "residue_numbering_policy": tables.residue_numbering_policy,
      "schema_version": geometry_bundle.SCHEMA_VERSION,
      "system_id": system_id,
      "thresholds": {
          "confidence_min_residues": MIN_CONFIDENCE_CORRELATION_RESIDUES,
          "ensemble_like_median_pairwise_backbone_rmsd": 0.5,
          "ensemble_like_p95_local_residue_variance": 0.25,
          "near_duplicate_median_pairwise_backbone_rmsd": 0.1,
          "near_duplicate_p95_local_residue_variance": 0.05,
          "shared_backbone_atoms_min": model_diversity.MIN_SHARED_BACKBONE_ATOMS,
          "shared_residues_min": model_diversity.MIN_SHARED_RESIDUES,
      },
  }


def _confidence_linkage_status(
    extraction_metrics: Sequence[Mapping[str, Any]],
) -> str:
  if all(row.get("confidence_atom_count_matches") for row in extraction_metrics):
    return "linked"
  return "mismatch"


def _mean_plddt_by_residue(
    atom_rows: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, str, str, str, str], float]:
  rank_1 = sorted({str(row["model_rank"]) for row in atom_rows})
  if not rank_1:
    return {}
  primary_rank = rank_1[0]
  grouped: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
  for row in atom_rows:
    if row["model_rank"] != primary_rank:
      continue
    plddt = row.get("plddt")
    if plddt is None or plddt == "":
      continue
    grouped[
        (
            str(row["system_id"]),
            str(row["chain_id"]),
            str(row["residue_number"]),
            str(row["insertion_code"]),
            str(row.get("residue_name", "")),
        )
    ].append(float(plddt))
  return {key: float(np.mean(values)) for key, values in grouped.items()}


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
  try:
    from scipy import stats

    result = stats.spearmanr(x, y)
    return float(result.statistic)
  except Exception:  # pylint: disable=broad-exception-caught
    x_rank = _rankdata(x)
    y_rank = _rankdata(y)
    if np.std(x_rank) == 0 or np.std(y_rank) == 0:
      return float("nan")
    return float(np.corrcoef(x_rank, y_rank)[0, 1])


def _rankdata(values: np.ndarray) -> np.ndarray:
  order = np.argsort(values, kind="mergesort")
  ranks = np.empty(len(values), dtype=float)
  ranks[order] = np.arange(1, len(values) + 1, dtype=float)
  return ranks


def _first_system_id(rows: Sequence[Mapping[str, Any]]) -> str:
  return str(rows[0]["system_id"]) if rows else ""
