#!/usr/bin/env python3
"""Export a deterministic AF3 geometry bundle from one AF3 job output."""

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
  parser.add_argument("--af3-output-dir", required=True, type=Path)
  parser.add_argument("--output-dir", required=True, type=Path)
  parser.add_argument("--system-id", required=True)
  parser.add_argument("--top-ranked-only", action="store_true")
  parser.add_argument("--emit-plots", action="store_true")
  args = parser.parse_args(argv)
  output_audit.export_geometry_bundle(
      af3_output_dir=args.af3_output_dir,
      output_dir=args.output_dir,
      system_id=args.system_id,
      top_ranked_only=args.top_ranked_only,
      emit_plots=args.emit_plots,
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
