from pathlib import Path

from tools import audit_af3_output_geometry_readiness
import nmr_test_utils


def test_audit_cli_writes_readiness_report(tmp_path: Path):
  nmr_test_utils.write_af3_job(tmp_path / "outputs", job_name="job_a", n_residues=6)

  exit_code = audit_af3_output_geometry_readiness.main(
      [
          "--af3-output-root",
          str(tmp_path / "outputs"),
          "--output-root",
          str(tmp_path / "audit"),
          "--emit-plots",
      ]
  )

  assert exit_code == 0
  assert (tmp_path / "audit" / "readiness.json").exists()
  assert (tmp_path / "audit" / "REPORT.md").exists()
  assert (tmp_path / "audit" / "figures" / "blocked_reason_counts.png").exists()
