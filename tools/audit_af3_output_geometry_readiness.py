#!/usr/bin/env python3
"""Audit AF3 output directories for mapping-grade geometry readiness."""

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
  parser.add_argument("--af3-output-root", required=True, type=Path)
  parser.add_argument("--output-root", required=True, type=Path)
  parser.add_argument("--emit-plots", action="store_true")
  args = parser.parse_args(argv)
  output_audit.audit_output_root(
      af3_output_root=args.af3_output_root,
      output_root=args.output_root,
      emit_plots=args.emit_plots,
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
