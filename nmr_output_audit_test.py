from pathlib import Path

from alphafold3.nmr import geometry_metrics
from alphafold3.nmr import output_audit
import nmr_test_utils


def test_audit_output_root_writes_deterministic_tables(tmp_path: Path):
  nmr_test_utils.write_af3_job(tmp_path / "outputs", job_name="job_a", n_residues=6)
  nmr_test_utils.write_af3_job(
      tmp_path / "outputs", job_name="job_b", n_residues=6, confidence_mismatch=True
  )

  summary = output_audit.audit_output_root(
      af3_output_root=tmp_path / "outputs",
      output_root=tmp_path / "audit",
      emit_plots=True,
  )

  assert summary["job_count"] == 2
  assert (tmp_path / "audit" / "bundle_readiness_summary.csv").exists()
  assert (tmp_path / "audit" / "readiness.json").exists()
  assert (tmp_path / "audit" / "REPORT.md").exists()
  assert (tmp_path / "audit" / "figures" / "readiness_verdict_counts.png").exists()


def test_compressed_cif_and_confidence_are_supported(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(
      tmp_path, n_residues=6, compressed=True
  )

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.manifest["readiness_verdict"] == "mapping_grade"
  assert (tmp_path / "bundle" / "af3_models" / "rank_001.cif").exists()


def test_confidence_agreement_requires_minimum_residue_count(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6, perturb=True)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert (
      result.confidence_rows[0]["confidence_geometry_agreement"]
      == geometry_metrics.CONFIDENCE_NOT_APPLICABLE
  )


def test_confidence_agreement_is_aligned_when_plddt_tracks_variance(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=35, perturb=True)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.confidence_rows[0]["n_residues"] >= 30
  assert (
      result.confidence_rows[0]["confidence_geometry_agreement"]
      == geometry_metrics.CONFIDENCE_ALIGNED
  )
