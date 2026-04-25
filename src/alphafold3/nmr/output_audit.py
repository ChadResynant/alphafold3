# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0.

"""AF3 output readiness audit and geometry bundle export orchestration."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
import csv
import dataclasses
import json
from pathlib import Path
from typing import Any

from alphafold3.nmr import geometry_bundle
from alphafold3.nmr import geometry_metrics
from alphafold3.nmr import model_diversity


@dataclasses.dataclass(frozen=True, slots=True)
class Af3JobOutput:
  job_dir: Path
  job_name: str
  source_models: tuple[geometry_bundle.SourceModel, ...]
  blocked_reasons: tuple[str, ...] = ()


@dataclasses.dataclass(frozen=True, slots=True)
class BundleExportResult:
  manifest: dict[str, Any]
  tables: geometry_bundle.BundleTables
  diversity: model_diversity.DiversityResult
  confidence_rows: list[dict[str, Any]]
  readiness_rows: list[dict[str, Any]]
  identity_rows: list[dict[str, Any]]
  completeness_rows: list[dict[str, Any]]


def discover_af3_jobs(path: Path) -> list[Af3JobOutput]:
  """Discovers AF3 job output directories under `path`."""
  path = path.resolve()
  if _looks_like_job_dir(path):
    return [discover_single_job(path)]
  jobs = [discover_single_job(child) for child in sorted(path.iterdir()) if child.is_dir() and _looks_like_job_dir(child)]
  return jobs


def discover_single_job(job_dir: Path) -> Af3JobOutput:
  job_dir = job_dir.resolve()
  ranking_files = sorted(job_dir.glob("*_ranking_scores.csv"))
  blocked: list[str] = []
  if ranking_files:
    ranking_file = ranking_files[0]
    job_name = ranking_file.name.removesuffix("_ranking_scores.csv")
    source_models = _models_from_ranking(job_dir, job_name, ranking_file)
  else:
    blocked.append(geometry_bundle.BLOCKED_MISSING_RANKED_MODELS)
    model_path = _find_one(job_dir, "*_model.cif") or _find_one(job_dir, "*_model.cif.zst")
    if model_path is None:
      source_models = []
      job_name = job_dir.name
    else:
      job_name = model_path.name.removesuffix("_model.cif").removesuffix("_model.cif.zst")
      source_models = [
          geometry_bundle.SourceModel(
              rank_index=1,
              model_rank="rank_001",
              seed="",
              sample="",
              ranking_score=None,
              model_path=model_path,
              confidence_path=_find_one(job_dir, f"{job_name}_confidences.json")
              or _find_one(job_dir, f"{job_name}_confidences.json.zst"),
              summary_confidence_path=_find_one(
                  job_dir, f"{job_name}_summary_confidences.json"
              )
              or _find_one(job_dir, f"{job_name}_summary_confidences.json.zst"),
          )
      ]
  missing_model_count = sum(1 for model in source_models if not model.model_path.exists())
  if missing_model_count:
    blocked.append(geometry_bundle.BLOCKED_MISSING_RANKED_MODELS)
  return Af3JobOutput(
      job_dir=job_dir,
      job_name=job_name,
      source_models=tuple(source_models),
      blocked_reasons=tuple(sorted(set(blocked))),
  )


def export_geometry_bundle(
    *,
    af3_output_dir: Path,
    output_dir: Path,
    system_id: str,
    top_ranked_only: bool = False,
    emit_plots: bool = False,
) -> BundleExportResult:
  """Exports a deterministic AF3 geometry bundle for one job output."""
  job = discover_single_job(af3_output_dir)
  source_models = list(job.source_models)
  if top_ranked_only:
    source_models = source_models[:1]
  tables = geometry_bundle.build_bundle_tables(
      system_id=system_id,
      source_models=source_models,
  )
  if job.blocked_reasons:
    combined_blockers = tuple(sorted(set(tables.blocked_reasons) | set(job.blocked_reasons)))
    tables = dataclasses.replace(
        tables,
        blocked_reasons=combined_blockers,
        readiness_verdict=geometry_bundle.readiness_verdict(combined_blockers),
        identity_grade=geometry_bundle.identity_grade(
            combined_blockers, tables.residue_numbering_policy
        ),
        geometry_grade=geometry_bundle.geometry_grade(combined_blockers),
    )
  diversity = model_diversity.compute_model_diversity(tables.atom_rows)
  confidence_rows = geometry_metrics.confidence_geometry_rows(
      tables.atom_rows, diversity.residue_variance_rows
  )
  readiness_rows = geometry_metrics.bundle_readiness_summary(
      system_id=system_id,
      tables=tables,
      diversity=diversity,
      confidence_rows=confidence_rows,
  )
  identity_rows = geometry_metrics.identity_metrics_rows(tables.extraction_metrics)
  completeness_rows = geometry_metrics.geometry_completeness_rows(tables.extraction_metrics)
  manifest = geometry_bundle.export_bundle_artifacts(
      bundle_dir=output_dir,
      system_id=system_id,
      source_output_dir=af3_output_dir,
      source_models=source_models,
      tables=tables,
      diversity_rows=diversity.summary_rows,
      local_class_rows=diversity.local_class_rows,
      pairwise_rows=diversity.pairwise_rows,
      residue_variance_rows=diversity.residue_variance_rows,
      confidence_rows=confidence_rows,
  )
  manifest["diversity_grade"] = diversity.diversity_grade
  manifest["confidence_geometry_agreement"] = (
      confidence_rows[0]["confidence_geometry_agreement"]
      if confidence_rows
      else geometry_metrics.CONFIDENCE_NOT_APPLICABLE
  )
  geometry_bundle.stable_json_dump(manifest, output_dir / "bundle_manifest.json")
  _write_analysis_tables(
      output_dir=output_dir,
      readiness_rows=readiness_rows,
      identity_rows=identity_rows,
      completeness_rows=completeness_rows,
      metrics=geometry_metrics.metrics_json(
          system_id=system_id,
          tables=tables,
          diversity=diversity,
          confidence_rows=confidence_rows,
      ),
  )
  figure_paths: list[str] = []
  if emit_plots:
    from alphafold3.nmr import geometry_figures

    figure_paths = geometry_figures.write_figures(
        output_dir=output_dir / "figures",
        readiness_rows=readiness_rows,
        blocked_reasons=tables.blocked_reasons,
        pairwise_rows=diversity.pairwise_rows,
        residue_variance_rows=diversity.residue_variance_rows,
        confidence_rows=confidence_rows,
        local_class_rows=diversity.local_class_rows,
    )
    figure_paths = [str(Path("figures") / Path(path).name) for path in figure_paths]
  geometry_bundle.write_report(
      path=output_dir / "REPORT.md",
      title="AF3 Geometry Bundle Readiness",
      manifest=manifest,
      metrics=geometry_metrics.metrics_json(
          system_id=system_id,
          tables=tables,
          diversity=diversity,
          confidence_rows=confidence_rows,
      ),
      figures=figure_paths,
  )
  return BundleExportResult(
      manifest=manifest,
      tables=tables,
      diversity=diversity,
      confidence_rows=confidence_rows,
      readiness_rows=readiness_rows,
      identity_rows=identity_rows,
      completeness_rows=completeness_rows,
  )


def audit_output_root(
    *,
    af3_output_root: Path,
    output_root: Path,
    emit_plots: bool = False,
) -> dict[str, Any]:
  """Audits every AF3 job output under a root and writes combined reports."""
  jobs = discover_af3_jobs(af3_output_root)
  output_root.mkdir(parents=True, exist_ok=True)
  readiness_rows: list[dict[str, Any]] = []
  identity_rows: list[dict[str, Any]] = []
  completeness_rows: list[dict[str, Any]] = []
  diversity_rows: list[dict[str, Any]] = []
  pairwise_rows: list[dict[str, Any]] = []
  residue_variance_rows: list[dict[str, Any]] = []
  confidence_rows: list[dict[str, Any]] = []
  local_class_rows: list[dict[str, Any]] = []
  for job in jobs:
    result = _evaluate_job(job)
    readiness_rows.extend(result["readiness_rows"])
    identity_rows.extend(result["identity_rows"])
    completeness_rows.extend(result["completeness_rows"])
    diversity_rows.extend(result["diversity_rows"])
    pairwise_rows.extend(result["pairwise_rows"])
    residue_variance_rows.extend(result["residue_variance_rows"])
    confidence_rows.extend(result["confidence_rows"])
    local_class_rows.extend(result["local_class_rows"])

  geometry_bundle.write_csv(output_root / "bundle_readiness_summary.csv", readiness_rows)
  geometry_bundle.write_csv(output_root / "identity_metrics.csv", identity_rows)
  geometry_bundle.write_csv(output_root / "geometry_completeness.csv", completeness_rows)
  geometry_bundle.write_csv(output_root / "geometry_diversity_summary.csv", diversity_rows)
  geometry_bundle.write_csv(output_root / "pairwise_model_rmsd.csv", pairwise_rows)
  geometry_bundle.write_csv(output_root / "per_residue_variance.csv", residue_variance_rows)
  geometry_bundle.write_csv(output_root / "confidence_geometry_correlation.csv", confidence_rows)
  geometry_bundle.write_csv(output_root / "local_class_diversity.csv", local_class_rows)

  summary = {
      "artifact_class": "af3_output_geometry_readiness_audit",
      "schema_version": geometry_bundle.SCHEMA_VERSION,
      "af3_output_root": str(af3_output_root),
      "job_count": len(jobs),
      "readiness_verdict_counts": dict(
          sorted(Counter(row["readiness_verdict"] for row in readiness_rows).items())
      ),
      "blocked_reason_counts": _blocked_reason_counts(readiness_rows),
      "promotion_allowed": False,
      "training_surface_allowed": False,
      "contains_nmr_residuals": False,
      "contains_reference_structure": False,
      "mapping_inference_allowed": False,
  }
  geometry_bundle.stable_json_dump(summary, output_root / "readiness.json")
  geometry_bundle.stable_json_dump(summary, output_root / "metrics.json")

  figure_paths: list[str] = []
  if emit_plots:
    from alphafold3.nmr import geometry_figures

    figure_paths = geometry_figures.write_figures(
        output_dir=output_root / "figures",
        readiness_rows=readiness_rows,
        blocked_reasons=[
            reason
            for row in readiness_rows
            for reason in str(row.get("blocked_reasons", "")).split(";")
            if reason
        ],
        pairwise_rows=pairwise_rows,
        residue_variance_rows=residue_variance_rows,
        confidence_rows=confidence_rows,
        local_class_rows=local_class_rows,
    )
    figure_paths = [str(Path("figures") / Path(path).name) for path in figure_paths]
  _write_audit_report(output_root / "REPORT.md", summary, figure_paths)
  return summary


def _evaluate_job(job: Af3JobOutput) -> dict[str, Any]:
  system_id = job.job_name
  tables = geometry_bundle.build_bundle_tables(
      system_id=system_id,
      source_models=job.source_models,
  )
  if job.blocked_reasons:
    combined_blockers = tuple(sorted(set(tables.blocked_reasons) | set(job.blocked_reasons)))
    tables = dataclasses.replace(
        tables,
        blocked_reasons=combined_blockers,
        readiness_verdict=geometry_bundle.readiness_verdict(combined_blockers),
        identity_grade=geometry_bundle.identity_grade(
            combined_blockers, tables.residue_numbering_policy
        ),
        geometry_grade=geometry_bundle.geometry_grade(combined_blockers),
    )
  diversity = model_diversity.compute_model_diversity(tables.atom_rows)
  confidence_rows = geometry_metrics.confidence_geometry_rows(
      tables.atom_rows, diversity.residue_variance_rows
  )
  return {
      "readiness_rows": geometry_metrics.bundle_readiness_summary(
          system_id=system_id,
          tables=tables,
          diversity=diversity,
          confidence_rows=confidence_rows,
      ),
      "identity_rows": geometry_metrics.identity_metrics_rows(tables.extraction_metrics),
      "completeness_rows": geometry_metrics.geometry_completeness_rows(
          tables.extraction_metrics
      ),
      "diversity_rows": diversity.summary_rows,
      "pairwise_rows": diversity.pairwise_rows,
      "residue_variance_rows": diversity.residue_variance_rows,
      "confidence_rows": confidence_rows,
      "local_class_rows": diversity.local_class_rows,
  }


def _write_analysis_tables(
    *,
    output_dir: Path,
    readiness_rows: Sequence[Mapping[str, Any]],
    identity_rows: Sequence[Mapping[str, Any]],
    completeness_rows: Sequence[Mapping[str, Any]],
    metrics: Mapping[str, Any],
) -> None:
  geometry_bundle.write_csv(output_dir / "bundle_readiness_summary.csv", readiness_rows)
  geometry_bundle.write_csv(output_dir / "identity_metrics.csv", identity_rows)
  geometry_bundle.write_csv(output_dir / "geometry_completeness.csv", completeness_rows)
  geometry_bundle.stable_json_dump(metrics, output_dir / "metrics.json")


def _models_from_ranking(
    job_dir: Path, job_name: str, ranking_file: Path
) -> list[geometry_bundle.SourceModel]:
  with ranking_file.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f))
  parsed = []
  for row in rows:
    score = _float_or_none(row.get("ranking_score"))
    parsed.append(
        {
            "seed": str(row.get("seed", "")),
            "sample": str(row.get("sample", "")),
            "ranking_score": score,
        }
    )
  parsed.sort(
      key=lambda row: (
          -(row["ranking_score"] if row["ranking_score"] is not None else float("-inf")),
          row["seed"],
          row["sample"],
      )
  )
  source_models = []
  for rank_index, row in enumerate(parsed, start=1):
    model_rank = f"rank_{rank_index:03d}"
    sample_dir = job_dir / f"seed-{row['seed']}_sample-{row['sample']}"
    sample_prefix = f"{job_name}_seed-{row['seed']}_sample-{row['sample']}"
    source_models.append(
        geometry_bundle.SourceModel(
            rank_index=rank_index,
            model_rank=model_rank,
            seed=row["seed"],
            sample=row["sample"],
            ranking_score=row["ranking_score"],
            model_path=_existing_or_default(
                sample_dir / f"{sample_prefix}_model.cif",
                sample_dir / f"{sample_prefix}_model.cif.zst",
            ),
            confidence_path=_find_existing(
                sample_dir / f"{sample_prefix}_confidences.json",
                sample_dir / f"{sample_prefix}_confidences.json.zst",
            ),
            summary_confidence_path=_find_existing(
                sample_dir / f"{sample_prefix}_summary_confidences.json",
                sample_dir / f"{sample_prefix}_summary_confidences.json.zst",
            ),
        )
    )
  return source_models


def _looks_like_job_dir(path: Path) -> bool:
  return bool(
      list(path.glob("*_ranking_scores.csv"))
      or list(path.glob("*_model.cif"))
      or list(path.glob("*_model.cif.zst"))
  )


def _find_one(path: Path, pattern: str) -> Path | None:
  matches = sorted(path.glob(pattern))
  return matches[0] if matches else None


def _find_existing(*paths: Path) -> Path | None:
  for path in paths:
    if path.exists():
      return path
  return None


def _existing_or_default(*paths: Path) -> Path:
  found = _find_existing(*paths)
  return found if found is not None else paths[0]


def _float_or_none(value: Any) -> float | None:
  try:
    if value in ("", None):
      return None
    return float(value)
  except ValueError:
    return None


def _blocked_reason_counts(rows: Sequence[Mapping[str, Any]]) -> dict[str, int]:
  counts: Counter[str] = Counter()
  for row in rows:
    for reason in str(row.get("blocked_reasons", "")).split(";"):
      if reason:
        counts[reason] += 1
  return dict(sorted(counts.items()))


def _write_audit_report(
    path: Path, summary: Mapping[str, Any], figures: Sequence[str]
) -> None:
  lines = [
      "# AF3 Output Geometry Readiness",
      "",
      f"- AF3 output root: `{summary['af3_output_root']}`",
      f"- Job count: `{summary['job_count']}`",
      "- `promotion_allowed`: false",
      "- `training_surface_allowed`: false",
      "- `contains_nmr_residuals`: false",
      "- `contains_reference_structure`: false",
      "- `mapping_inference_allowed`: false",
      "",
      "## Readiness Verdict Counts",
      "",
  ]
  for verdict, count in summary["readiness_verdict_counts"].items():
    lines.append(f"- `{verdict}`: `{count}`")
  lines.extend(["", "## Blocked Reason Counts", ""])
  if summary["blocked_reason_counts"]:
    for reason, count in summary["blocked_reason_counts"].items():
      lines.append(f"- `{reason}`: `{count}`")
  else:
    lines.append("- none")
  if figures:
    lines.extend(["", "## Figures", ""])
    for figure in figures:
      lines.append(f"![{Path(figure).stem}]({figure})")
  path.write_text("\n".join(lines) + "\n", encoding="utf-8")
