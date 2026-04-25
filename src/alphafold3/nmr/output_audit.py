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


BLOCKED_NO_REAL_AF3_OUTPUTS = "BLOCKED_NO_REAL_AF3_OUTPUTS"


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


def discover_real_af3_jobs(
    search_roots: Sequence[Path],
    *,
    exclude_paths: Sequence[Path] = (),
) -> list[Af3JobOutput]:
  """Recursively locates candidate real AF3 job output directories.

  This intentionally searches only for AF3 output-layout artifacts. It does not
  infer or repair output directories from filenames beyond identifying the AF3
  job root that owns a ranking CSV or model CIF.
  """
  excluded = [path.resolve() for path in exclude_paths]
  job_dirs: set[Path] = set()
  for root in search_roots:
    root = root.resolve()
    if not root.exists() or _is_excluded(root, excluded):
      continue
    if root.is_file():
      continue
    if _looks_like_job_dir(root):
      job_dirs.add(root)
      continue
    patterns = ("*_ranking_scores.csv", "*_model.cif", "*_model.cif.zst")
    for pattern in patterns:
      for artifact in root.rglob(pattern):
        if _is_excluded(artifact, excluded):
          continue
        job_dir = _job_dir_from_artifact(artifact)
        if job_dir is not None and not _is_excluded(job_dir, excluded):
          job_dirs.add(job_dir)
  return [discover_single_job(path) for path in sorted(job_dirs)]


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


def survey_real_af3_outputs(
    *,
    search_roots: Sequence[Path],
    output_root: Path,
    exclude_paths: Sequence[Path] = (),
    emit_plots: bool = False,
) -> dict[str, Any]:
  """Runs the Phase D real AF3 geometry-readiness cohort survey."""
  output_root.mkdir(parents=True, exist_ok=True)
  excludes = tuple([output_root, *exclude_paths])
  jobs = discover_real_af3_jobs(search_roots, exclude_paths=excludes)
  if not jobs:
    summary = _write_blocked_real_survey(
        search_roots=search_roots,
        output_root=output_root,
        exclude_paths=excludes,
    )
    return summary

  metrics_rows: list[dict[str, Any]] = []
  for job in jobs:
    bundle_dir = output_root / "bundles" / job.job_name
    if job.source_models and all(model.model_path.exists() for model in job.source_models):
      result = export_geometry_bundle(
          af3_output_dir=job.job_dir,
          output_dir=bundle_dir,
          system_id=job.job_name,
          emit_plots=emit_plots,
      )
      metrics_rows.extend(
          _cohort_rows_from_result(job=job, result=result)
      )
    else:
      result = _evaluate_job(job)
      metrics_rows.extend(_cohort_rows_from_evaluation(job=job, result=result))

  blocked_rows = _count_rows_from_semicolon_column(
      metrics_rows,
      source_column="blocked_reasons",
      output_name="blocked_reason",
  )
  diversity_rows = _count_rows(
      metrics_rows,
      source_column="diversity_grade",
      output_name="diversity_grade",
  )
  confidence_rows = _count_rows(
      metrics_rows,
      source_column="confidence_geometry_agreement",
      output_name="confidence_geometry_agreement",
  )
  geometry_rows = _count_rows(
      metrics_rows,
      source_column="geometry_grade",
      output_name="geometry_grade",
  )
  readiness_rows = _count_rows(
      metrics_rows,
      source_column="readiness_verdict",
      output_name="readiness_verdict",
  )
  geometry_bundle.write_csv(output_root / "metrics_cohort.csv", metrics_rows)
  geometry_bundle.write_csv(output_root / "blocked_reasons.csv", blocked_rows)
  geometry_bundle.write_csv(
      output_root / "diversity_grade_summary.csv", diversity_rows
  )
  geometry_bundle.write_csv(
      output_root / "confidence_geometry_summary.csv", confidence_rows
  )
  summary = {
      "artifact_class": "af3_real_geometry_readiness_survey",
      "blocked_reason_counts": {
          row["blocked_reason"]: row["count"] for row in blocked_rows
      },
      "confidence_geometry_agreement_counts": {
          row["confidence_geometry_agreement"]: row["count"]
          for row in confidence_rows
      },
      "contains_nmr_residuals": False,
      "contains_reference_structure": False,
      "diversity_grade_counts": {
          row["diversity_grade"]: row["count"] for row in diversity_rows
      },
      "geometry_grade_counts": {
          row["geometry_grade"]: row["count"] for row in geometry_rows
      },
      "job_count": len(metrics_rows),
      "mapping_inference_allowed": False,
      "promotion_allowed": False,
      "readiness_verdict_counts": {
          row["readiness_verdict"]: row["count"] for row in readiness_rows
      },
      "schema_version": geometry_bundle.SCHEMA_VERSION,
      "search_roots": [str(path) for path in search_roots],
      "status": "ok",
      "training_surface_allowed": False,
  }
  geometry_bundle.stable_json_dump(summary, output_root / "metrics.json")
  _write_real_survey_report(
      output_root / "REPORT.md",
      summary=summary,
      blocked_rows=blocked_rows,
      diversity_rows=diversity_rows,
      confidence_rows=confidence_rows,
      geometry_rows=geometry_rows,
      readiness_rows=readiness_rows,
  )
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


def _cohort_rows_from_result(
    *,
    job: Af3JobOutput,
    result: BundleExportResult,
) -> list[dict[str, Any]]:
  readiness = result.readiness_rows[0] if result.readiness_rows else {}
  confidence = result.confidence_rows[0] if result.confidence_rows else {}
  diversity = result.diversity.summary_rows[0] if result.diversity.summary_rows else {}
  return [
      {
          "system_id": job.job_name,
          "source_output_dir": str(job.job_dir),
          "readiness_verdict": result.manifest.get("readiness_verdict", ""),
          "identity_grade": result.manifest.get("identity_grade", ""),
          "geometry_grade": result.manifest.get("geometry_grade", ""),
          "diversity_grade": result.diversity.diversity_grade,
          "confidence_geometry_agreement": confidence.get(
              "confidence_geometry_agreement",
              geometry_metrics.CONFIDENCE_NOT_APPLICABLE,
          ),
          "blocked_reasons": ";".join(result.manifest.get("blocked_reasons", [])),
          "n_models": readiness.get("n_models", ""),
          "n_chains": readiness.get("n_chains", ""),
          "n_residues": readiness.get("n_residues", ""),
          "n_atoms": readiness.get("n_atoms", ""),
          "median_pairwise_backbone_rmsd": diversity.get(
              "median_pairwise_backbone_rmsd", ""
          ),
          "p95_per_residue_variance": diversity.get(
              "p95_per_residue_variance", ""
          ),
          "confidence_spearman_plddt_vs_variance": confidence.get(
              "spearman_plddt_vs_variance", ""
          ),
          "atom_table_hash": result.manifest.get("atom_table_hash", ""),
          "bundle_dir": str(Path("bundles") / job.job_name),
      }
  ]


def _cohort_rows_from_evaluation(
    *, job: Af3JobOutput, result: Mapping[str, Any]
) -> list[dict[str, Any]]:
  readiness = result["readiness_rows"][0] if result["readiness_rows"] else {}
  confidence = result["confidence_rows"][0] if result["confidence_rows"] else {}
  diversity = result["diversity_rows"][0] if result["diversity_rows"] else {}
  return [
      {
          "system_id": job.job_name,
          "source_output_dir": str(job.job_dir),
          "readiness_verdict": readiness.get(
              "readiness_verdict", geometry_bundle.READINESS_NOT_MAPPING_GRADE
          ),
          "identity_grade": readiness.get("identity_grade", ""),
          "geometry_grade": readiness.get("geometry_grade", ""),
          "diversity_grade": diversity.get(
              "diversity_grade", model_diversity.DIVERSITY_NOT_APPLICABLE
          ),
          "confidence_geometry_agreement": confidence.get(
              "confidence_geometry_agreement",
              geometry_metrics.CONFIDENCE_NOT_APPLICABLE,
          ),
          "blocked_reasons": readiness.get("blocked_reasons", ""),
          "n_models": readiness.get("n_models", ""),
          "n_chains": readiness.get("n_chains", ""),
          "n_residues": readiness.get("n_residues", ""),
          "n_atoms": readiness.get("n_atoms", ""),
          "median_pairwise_backbone_rmsd": diversity.get(
              "median_pairwise_backbone_rmsd", ""
          ),
          "p95_per_residue_variance": diversity.get(
              "p95_per_residue_variance", ""
          ),
          "confidence_spearman_plddt_vs_variance": confidence.get(
              "spearman_plddt_vs_variance", ""
          ),
          "atom_table_hash": "",
          "bundle_dir": "",
      }
  ]


def _write_blocked_real_survey(
    *,
    search_roots: Sequence[Path],
    output_root: Path,
    exclude_paths: Sequence[Path],
) -> dict[str, Any]:
  metrics_rows = [
      {
          "status": BLOCKED_NO_REAL_AF3_OUTPUTS,
          "job_count": 0,
          "search_roots": ";".join(str(path) for path in search_roots),
          "excluded_paths": ";".join(str(path) for path in exclude_paths),
      }
  ]
  blocked_rows = [{"blocked_reason": BLOCKED_NO_REAL_AF3_OUTPUTS, "count": 1}]
  diversity_rows = [
      {"diversity_grade": model_diversity.DIVERSITY_NOT_APPLICABLE, "count": 0}
  ]
  confidence_rows = [
      {
          "confidence_geometry_agreement": geometry_metrics.CONFIDENCE_NOT_APPLICABLE,
          "count": 0,
      }
  ]
  geometry_bundle.write_csv(output_root / "metrics_cohort.csv", metrics_rows)
  geometry_bundle.write_csv(output_root / "blocked_reasons.csv", blocked_rows)
  geometry_bundle.write_csv(
      output_root / "diversity_grade_summary.csv", diversity_rows
  )
  geometry_bundle.write_csv(
      output_root / "confidence_geometry_summary.csv", confidence_rows
  )
  summary = {
      "artifact_class": "af3_real_geometry_readiness_survey",
      "blocked_reason": BLOCKED_NO_REAL_AF3_OUTPUTS,
      "blocked_reason_counts": {BLOCKED_NO_REAL_AF3_OUTPUTS: 1},
      "confidence_geometry_agreement_counts": {},
      "contains_nmr_residuals": False,
      "contains_reference_structure": False,
      "diversity_grade_counts": {},
      "geometry_grade_counts": {},
      "job_count": 0,
      "mapping_inference_allowed": False,
      "promotion_allowed": False,
      "readiness_verdict_counts": {},
      "schema_version": geometry_bundle.SCHEMA_VERSION,
      "search_roots": [str(path) for path in search_roots],
      "status": BLOCKED_NO_REAL_AF3_OUTPUTS,
      "training_surface_allowed": False,
  }
  geometry_bundle.stable_json_dump(summary, output_root / "metrics.json")
  _write_blocked_real_survey_report(output_root / "REPORT.md", summary)
  return summary


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


def _is_excluded(path: Path, excluded: Sequence[Path]) -> bool:
  resolved = path.resolve()
  for excluded_path in excluded:
    try:
      resolved.relative_to(excluded_path)
      return True
    except ValueError:
      continue
  return any(part in {".git", "__pycache__"} for part in resolved.parts)


def _job_dir_from_artifact(path: Path) -> Path | None:
  if path.name.endswith("_ranking_scores.csv"):
    return path.parent
  if path.name.endswith("_model.cif") or path.name.endswith("_model.cif.zst"):
    parent = path.parent
    if parent.name.startswith("seed-") and "_sample-" in parent.name:
      return parent.parent
    return parent
  return None


def _count_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    source_column: str,
    output_name: str,
) -> list[dict[str, Any]]:
  counts = Counter(str(row.get(source_column, "")) for row in rows)
  return [
      {output_name: value, "count": count}
      for value, count in sorted(counts.items())
      if value
  ]


def _count_rows_from_semicolon_column(
    rows: Sequence[Mapping[str, Any]],
    *,
    source_column: str,
    output_name: str,
) -> list[dict[str, Any]]:
  counts: Counter[str] = Counter()
  for row in rows:
    values = [value for value in str(row.get(source_column, "")).split(";") if value]
    if not values:
      values = ["none"]
    counts.update(values)
  return [
      {output_name: value, "count": count}
      for value, count in sorted(counts.items())
  ]


def _write_real_survey_report(
    path: Path,
    *,
    summary: Mapping[str, Any],
    blocked_rows: Sequence[Mapping[str, Any]],
    diversity_rows: Sequence[Mapping[str, Any]],
    confidence_rows: Sequence[Mapping[str, Any]],
    geometry_rows: Sequence[Mapping[str, Any]],
    readiness_rows: Sequence[Mapping[str, Any]],
) -> None:
  lines = [
      "# Phase D Real AF3 Geometry Readiness Survey",
      "",
      f"- Status: `{summary['status']}`",
      f"- Job count: `{summary['job_count']}`",
      "- `promotion_allowed`: false",
      "- `training_surface_allowed`: false",
      "- `contains_nmr_residuals`: false",
      "- `contains_reference_structure`: false",
      "- `mapping_inference_allowed`: false",
      "",
      "## Mapping Readiness",
      "",
  ]
  lines.extend(_markdown_count_lines(readiness_rows, "readiness_verdict"))
  lines.extend(["", "## Geometry Grades", ""])
  lines.extend(_markdown_count_lines(geometry_rows, "geometry_grade"))
  lines.extend(["", "## Diversity Grades", ""])
  lines.extend(_markdown_count_lines(diversity_rows, "diversity_grade"))
  lines.extend(["", "## Confidence/Geometry Agreement", ""])
  lines.extend(
      _markdown_count_lines(confidence_rows, "confidence_geometry_agreement")
  )
  lines.extend(["", "## Blocked Reasons", ""])
  lines.extend(_markdown_count_lines(blocked_rows, "blocked_reason"))
  lines.extend(
      [
          "",
          "## Output Tables",
          "",
          "- `metrics_cohort.csv`",
          "- `blocked_reasons.csv`",
          "- `diversity_grade_summary.csv`",
          "- `confidence_geometry_summary.csv`",
      ]
  )
  path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_blocked_real_survey_report(
    path: Path, summary: Mapping[str, Any]
) -> None:
  lines = [
      "# Phase D Real AF3 Geometry Readiness Survey",
      "",
      f"- Status: `{BLOCKED_NO_REAL_AF3_OUTPUTS}`",
      "- Job count: `0`",
      "- `promotion_allowed`: false",
      "- `training_surface_allowed`: false",
      "- `contains_nmr_residuals`: false",
      "- `contains_reference_structure`: false",
      "- `mapping_inference_allowed`: false",
      "",
      "## Blocker",
      "",
      f"- `{BLOCKED_NO_REAL_AF3_OUTPUTS}`",
      "",
      "No real AF3 output directories were found under the configured search roots.",
      "",
      "## Search Roots",
      "",
  ]
  for root in summary["search_roots"]:
    lines.append(f"- `{root}`")
  path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _markdown_count_lines(
    rows: Sequence[Mapping[str, Any]], value_key: str
) -> list[str]:
  if not rows:
    return ["- none"]
  return [f"- `{row[value_key]}`: `{row['count']}`" for row in rows]


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
