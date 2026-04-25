from pathlib import Path

from alphafold3.nmr import model_diversity
from alphafold3.nmr import output_audit
import nmr_test_utils


def test_near_duplicate_ranked_models_are_classified(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6, perturb=False)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.diversity.diversity_grade == model_diversity.DIVERSITY_NEAR_DUPLICATE


def test_perturbed_ranked_models_are_ensemble_like(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=8, perturb=True)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.diversity.diversity_grade == model_diversity.DIVERSITY_ENSEMBLE_LIKE
  assert result.diversity.p95_per_residue_variance is not None


def test_insufficient_shared_data_is_not_applicable(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=4, perturb=True)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.diversity.diversity_grade == model_diversity.DIVERSITY_NOT_APPLICABLE


def test_amide_vector_requires_explicit_hydrogen(tmp_path: Path):
  no_h_job = nmr_test_utils.write_af3_job(tmp_path / "no_h", n_residues=6)
  with_h_job = nmr_test_utils.write_af3_job(
      tmp_path / "with_h", n_residues=6, hydrogens=True, perturb=True
  )

  no_h = output_audit.export_geometry_bundle(
      af3_output_dir=no_h_job,
      output_dir=tmp_path / "no_h_bundle",
      system_id="no_h",
  )
  with_h = output_audit.export_geometry_bundle(
      af3_output_dir=with_h_job,
      output_dir=tmp_path / "with_h_bundle",
      system_id="with_h",
  )

  no_h_amide = [
      row for row in no_h.diversity.local_class_rows if row["class"] == "amide_vector_angle"
  ][0]
  with_h_amide = [
      row for row in with_h.diversity.local_class_rows if row["class"] == "amide_vector_angle"
  ][0]
  assert no_h_amide["status"] == "not_applicable"
  assert with_h_amide["status"] == "ok"
