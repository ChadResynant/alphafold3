#!/usr/bin/env python3
"""Run the Phase D real AF3 geometry-readiness cohort survey."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
  sys.path.insert(0, str(SRC_ROOT))

from alphafold3.nmr import output_audit  # noqa: E402


def main(argv: list[str] | None = None) -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
      "--search-root",
      action="append",
      type=Path,
      default=[],
      help="Root to recursively search for real AF3 output directories.",
  )
  parser.add_argument(
      "--output-root",
      type=Path,
      default=Path("artifacts/nmr_geometry_readiness/phase_d_real_survey"),
  )
  parser.add_argument(
      "--exclude-path",
      action="append",
      type=Path,
      default=[],
      help="Path to exclude from recursive discovery.",
  )
  parser.add_argument("--emit-plots", action="store_true")
  args = parser.parse_args(argv)
  search_roots = args.search_root or [Path(".")]
  output_audit.survey_real_af3_outputs(
      search_roots=search_roots,
      output_root=args.output_root,
      exclude_paths=args.exclude_path,
      emit_plots=args.emit_plots,
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
