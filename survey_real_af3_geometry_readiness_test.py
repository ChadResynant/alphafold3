from pathlib import Path

from alphafold3.nmr import output_audit
from tools import survey_real_af3_geometry_readiness
import nmr_test_utils


def test_real_survey_blocks_when_no_real_outputs(tmp_path: Path):
  summary = output_audit.survey_real_af3_outputs(
      search_roots=[tmp_path / "empty"],
      output_root=tmp_path / "survey",
  )

  assert summary["status"] == output_audit.BLOCKED_NO_REAL_AF3_OUTPUTS
  assert (tmp_path / "survey" / "REPORT.md").exists()
  assert (tmp_path / "survey" / "metrics_cohort.csv").exists()
  assert (
      output_audit.BLOCKED_NO_REAL_AF3_OUTPUTS
      in (tmp_path / "survey" / "REPORT.md").read_text()
  )


def test_real_survey_writes_cohort_summaries_for_discovered_outputs(
    tmp_path: Path,
):
  nmr_test_utils.write_af3_job(
      tmp_path / "real_outputs",
      job_name="real_a",
      n_residues=35,
      n_models=3,
      perturb=True,
  )
  nmr_test_utils.write_af3_job(
      tmp_path / "real_outputs",
      job_name="real_b",
      n_residues=6,
      n_models=2,
      confidence_mismatch=True,
  )

  summary = output_audit.survey_real_af3_outputs(
      search_roots=[tmp_path / "real_outputs"],
      output_root=tmp_path / "survey",
  )

  assert summary["status"] == "ok"
  assert summary["job_count"] == 2
  assert summary["readiness_verdict_counts"]["mapping_grade"] == 1
  assert summary["readiness_verdict_counts"]["not_mapping_grade"] == 1
  assert (tmp_path / "survey" / "metrics_cohort.csv").exists()
  assert (tmp_path / "survey" / "blocked_reasons.csv").exists()
  assert (tmp_path / "survey" / "diversity_grade_summary.csv").exists()
  assert (tmp_path / "survey" / "confidence_geometry_summary.csv").exists()
  assert (tmp_path / "survey" / "bundles" / "real_a" / "REPORT.md").exists()


def test_real_survey_cli_writes_blocked_artifact(tmp_path: Path):
  exit_code = survey_real_af3_geometry_readiness.main(
      [
          "--search-root",
          str(tmp_path / "empty"),
          "--output-root",
          str(tmp_path / "survey"),
      ]
  )

  assert exit_code == 0
  assert (tmp_path / "survey" / "REPORT.md").exists()
  assert (
      output_audit.BLOCKED_NO_REAL_AF3_OUTPUTS
      in (tmp_path / "survey" / "metrics.json").read_text()
  )
