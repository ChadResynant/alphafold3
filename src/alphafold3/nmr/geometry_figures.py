# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0.

"""Deterministic PNG figures for AF3 geometry readiness reports."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


def write_figures(
    *,
    output_dir: Path,
    readiness_rows: Sequence[Mapping[str, Any]],
    blocked_reasons: Sequence[str],
    pairwise_rows: Sequence[Mapping[str, Any]],
    residue_variance_rows: Sequence[Mapping[str, Any]],
    confidence_rows: Sequence[Mapping[str, Any]],
    local_class_rows: Sequence[Mapping[str, Any]],
) -> list[str]:
  """Writes deterministic PNG figures and returns their paths."""
  plt = _matplotlib()
  output_dir.mkdir(parents=True, exist_ok=True)
  paths = [
      _bar(
          plt,
          output_dir / "readiness_verdict_counts.png",
          Counter(row.get("readiness_verdict", "") for row in readiness_rows),
          "Readiness verdict counts",
          "verdict",
      ),
      _bar(
          plt,
          output_dir / "blocked_reason_counts.png",
          Counter(reason for reason in blocked_reasons if reason),
          "Blocked reason counts",
          "reason",
      ),
      _pairwise_histogram(
          plt, output_dir / "pairwise_rmsd_histogram.png", pairwise_rows
      ),
      _rank_vs_rank1(
          plt, output_dir / "rank_vs_rank1_backbone_rmsd.png", pairwise_rows
      ),
      _residue_variance_trace(
          plt, output_dir / "per_residue_variance_trace.png", residue_variance_rows
      ),
      _confidence_scatter(
          plt, output_dir / "plddt_vs_coordinate_variance.png", confidence_rows
      ),
      _confidence_scatter(
          plt,
          output_dir / "confidence_binned_variance_distributions.png",
          confidence_rows,
      ),
      _sidechain_vs_backbone(
          plt,
          output_dir / "sidechain_vs_backbone_variance.png",
          residue_variance_rows,
      ),
      _local_class_bar(
          plt, output_dir / "local_class_diversity.png", local_class_rows
      ),
      _class_single_bar(
          plt,
          output_dir / "aromatic_ring_centroid_plane_variance.png",
          local_class_rows,
          ("aromatic_ring_centroid",),
      ),
      _class_single_bar(
          plt,
          output_dir / "methyl_carbon_variance.png",
          local_class_rows,
          ("methyl_carbon",),
      ),
      _class_single_bar(
          plt,
          output_dir / "amide_vector_variance.png",
          local_class_rows,
          ("amide_vector_angle",),
      ),
  ]
  return [str(path) for path in paths if path is not None]


def _matplotlib():
  import matplotlib

  matplotlib.use("Agg", force=True)
  import matplotlib.pyplot as plt

  plt.rcParams.update(
      {
          "axes.grid": True,
          "axes.titlesize": 11,
          "axes.labelsize": 9,
          "font.family": "DejaVu Sans",
          "figure.dpi": 120,
          "savefig.dpi": 120,
          "legend.fontsize": 8,
          "xtick.labelsize": 8,
          "ytick.labelsize": 8,
      }
  )
  return plt


def _bar(plt, path: Path, counts: Counter[str], title: str, xlabel: str) -> Path:
  labels = sorted(label for label in counts if label)
  values = [counts[label] for label in labels]
  if not labels:
    labels = ["none"]
    values = [0]
  fig, ax = plt.subplots(figsize=(6, 3.5), dpi=120)
  ax.bar(labels, values, color="#4c78a8")
  ax.set_title(title)
  ax.set_xlabel(xlabel)
  ax.set_ylabel("count")
  ax.tick_params(axis="x", labelrotation=30)
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _pairwise_histogram(plt, path: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
  values = [
      float(row["backbone_rmsd"])
      for row in rows
      if row.get("backbone_rmsd") not in ("", None)
  ]
  fig, ax = plt.subplots(figsize=(6, 3.5), dpi=120)
  ax.hist(values or [0.0], bins=min(20, max(1, len(values))), color="#59a14f")
  ax.set_title("Pairwise backbone RMSD")
  ax.set_xlabel("RMSD (A)")
  ax.set_ylabel("model pairs")
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _residue_variance_trace(
    plt, path: Path, rows: Sequence[Mapping[str, Any]]
) -> Path:
  sorted_rows = sorted(
      rows,
      key=lambda row: (
          str(row.get("chain_id", "")),
          _residue_number_sort_value(row.get("residue_number", "")),
          str(row.get("insertion_code", "")),
      ),
  )
  values = [
      float(row["coordinate_variance"])
      for row in sorted_rows
      if row.get("coordinate_variance") not in ("", None)
  ]
  fig, ax = plt.subplots(figsize=(7, 3.5), dpi=120)
  ax.plot(range(1, len(values) + 1), values, color="#f28e2b", linewidth=1.5)
  ax.set_title("Per-residue coordinate variance")
  ax.set_xlabel("ordered residue")
  ax.set_ylabel("variance (A^2)")
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _rank_vs_rank1(plt, path: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
  rank_rows = [
      row
      for row in rows
      if row.get("model_rank_a") == "rank_001"
      and row.get("backbone_rmsd") not in ("", None)
  ]
  rank_rows = sorted(rank_rows, key=lambda row: str(row.get("model_rank_b", "")))
  values = [float(row["backbone_rmsd"]) for row in rank_rows]
  labels = [str(row["model_rank_b"]) for row in rank_rows]
  fig, ax = plt.subplots(figsize=(6, 3.5), dpi=120)
  ax.plot(labels, values, marker="o", color="#4e79a7")
  ax.set_title("Rank versus rank-1 backbone RMSD")
  ax.set_xlabel("model rank")
  ax.set_ylabel("RMSD (A)")
  ax.tick_params(axis="x", labelrotation=30)
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _confidence_scatter(
    plt, path: Path, rows: Sequence[Mapping[str, Any]]
) -> Path:
  # The correlation table is intentionally compact; this figure records the
  # aggregate status rather than reconstructing raw per-residue pairs.
  labels = [str(row.get("confidence_geometry_agreement", "")) for row in rows]
  counts = Counter(label for label in labels if label)
  return _bar(plt, path, counts, "Confidence/geometry agreement", "agreement")


def _sidechain_vs_backbone(
    plt, path: Path, rows: Sequence[Mapping[str, Any]]
) -> Path:
  backbone = [
      float(row["backbone_variance"])
      for row in rows
      if row.get("backbone_variance") not in ("", None)
  ]
  sidechain = [
      float(row["sidechain_variance"])
      for row in rows
      if row.get("sidechain_variance") not in ("", None)
  ]
  labels = ["backbone", "sidechain"]
  values = [
      sum(backbone) / len(backbone) if backbone else 0.0,
      sum(sidechain) / len(sidechain) if sidechain else 0.0,
  ]
  fig, ax = plt.subplots(figsize=(5, 3.5), dpi=120)
  ax.bar(labels, values, color=["#4e79a7", "#f28e2b"])
  ax.set_title("Sidechain versus backbone variance")
  ax.set_ylabel("mean variance (A^2)")
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _local_class_bar(
    plt, path: Path, rows: Sequence[Mapping[str, Any]]
) -> Path:
  ok_rows = [row for row in rows if row.get("mean_variance") not in ("", None)]
  labels = [str(row.get("class", "")) for row in ok_rows]
  values = [float(row["mean_variance"]) for row in ok_rows]
  if not labels:
    labels = ["not_applicable"]
    values = [0.0]
  fig, ax = plt.subplots(figsize=(7, 3.5), dpi=120)
  ax.bar(labels, values, color="#b07aa1")
  ax.set_title("Local structural-class diversity")
  ax.set_xlabel("class")
  ax.set_ylabel("mean variance")
  ax.tick_params(axis="x", labelrotation=30)
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _class_single_bar(
    plt,
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    classes: tuple[str, ...],
) -> Path:
  selected = [row for row in rows if row.get("class") in classes]
  labels = [str(row.get("class", "")) for row in selected] or [classes[0]]
  values = [
      float(row["mean_variance"])
      if row.get("mean_variance") not in ("", None)
      else 0.0
      for row in selected
  ] or [0.0]
  fig, ax = plt.subplots(figsize=(5, 3.5), dpi=120)
  ax.bar(labels, values, color="#76b7b2")
  ax.set_title(classes[0].replace("_", " "))
  ax.set_ylabel("mean variance")
  ax.tick_params(axis="x", labelrotation=20)
  fig.tight_layout()
  _savefig(fig, path)
  plt.close(fig)
  return path


def _savefig(fig, path: Path) -> None:
  fig.savefig(
      path,
      format="png",
      metadata={
          "Software": "",
          "Creation Time": "",
          "Description": "",
      },
  )


def _residue_number_sort_value(value: Any) -> tuple[int, Any]:
  text = str(value)
  try:
    return (0, int(text))
  except ValueError:
    return (1, text)
