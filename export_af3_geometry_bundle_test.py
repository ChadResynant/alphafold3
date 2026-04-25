from pathlib import Path

from tools import export_af3_geometry_bundle
import nmr_test_utils


def test_export_cli_writes_bundle(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6)

  exit_code = export_af3_geometry_bundle.main(
      [
          "--af3-output-dir",
          str(job_dir),
          "--output-dir",
          str(tmp_path / "bundle"),
          "--system-id",
          "sys1",
          "--emit-plots",
      ]
  )

  assert exit_code == 0
  assert (tmp_path / "bundle" / "bundle_manifest.json").exists()
  assert (tmp_path / "bundle" / "atom_identity_table.parquet").exists()
  assert (tmp_path / "bundle" / "figures" / "pairwise_rmsd_histogram.png").exists()


def test_export_cli_top_ranked_only_makes_diversity_not_applicable(tmp_path: Path):
  job_dir = nmr_test_utils.write_af3_job(tmp_path, n_residues=6)

  export_af3_geometry_bundle.main(
      [
          "--af3-output-dir",
          str(job_dir),
          "--output-dir",
          str(tmp_path / "bundle"),
          "--system-id",
          "sys1",
          "--top-ranked-only",
      ]
  )

  manifest = (tmp_path / "bundle" / "bundle_manifest.json").read_text()
  assert '"model_count": 1' in manifest
  assert '"diversity_grade": "not_applicable"' in manifest
