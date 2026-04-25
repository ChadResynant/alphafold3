from pathlib import Path

from alphafold3.nmr import geometry_bundle
from alphafold3.nmr import output_audit
import nmr_test_utils


def test_geometry_bundle_extracts_canonical_identity_with_element(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6)

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.manifest["readiness_verdict"] == "mapping_grade"
  assert result.manifest["mapping_inference_allowed"] is False
  assert result.manifest["atom_table_hash"]
  first_atom = result.tables.atom_rows[0]
  assert first_atom["element"]
  assert first_atom["atom_identity_hash"] == geometry_bundle.atom_identity_hash(
      first_atom
  )
  assert (tmp_path / "bundle" / "atom_identity_table.parquet").exists()


def test_duplicate_canonical_atom_identity_fails_closed(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(
      tmp_path, n_residues=6, duplicate_atom=True
  )

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.manifest["readiness_verdict"] == "not_mapping_grade"
  assert (
      geometry_bundle.BLOCKED_DUPLICATE_ATOM_IDENTITIES
      in result.manifest["blocked_reasons"]
  )


def test_disallowed_altloc_fails_closed(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6, altloc="B")

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.manifest["readiness_verdict"] == "not_mapping_grade"
  assert geometry_bundle.BLOCKED_NONTRIVIAL_ALTLOC in result.manifest["blocked_reasons"]


def test_confidence_atom_count_mismatch_fails_closed(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(
      tmp_path, n_residues=6, confidence_mismatch=True
  )

  result = output_audit.export_geometry_bundle(
      af3_output_dir=job_dir,
      output_dir=tmp_path / "bundle",
      system_id="sys1",
  )

  assert result.manifest["readiness_verdict"] == "not_mapping_grade"
  assert (
      geometry_bundle.BLOCKED_CONFIDENCE_ATOM_COUNT_MISMATCH
      in result.manifest["blocked_reasons"]
  )
