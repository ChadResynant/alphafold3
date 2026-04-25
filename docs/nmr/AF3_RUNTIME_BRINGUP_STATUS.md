# AF3 Runtime Bring-Up Status

Phase F goal: produce one valid real AF3 output directory, then run the Phase C
geometry-readiness audit/export on that output.

Current status: `BLOCKED_AF3_RUNTIME_PREREQUISITES_ABSENT`.

See [AF3_RUNTIME_UNBLOCK_CHECKLIST.md](AF3_RUNTIME_UNBLOCK_CHECKLIST.md)
for the one-page “proper machine” checklist to get to one real golden-path run.

## What Was Checked

- Active Python is `/Users/chad/.pyenv/shims/python`, version `3.11.7`.
- `absl-py` was installed and now imports in the active Python.
- Active Python cannot import `lzma` because `_lzma` is missing.
- Homebrew Python 3.13 can import `lzma`, but lacks the AF3 dependency stack
  needed for this checkout.
- No AF3 model parameter directory was found in the bounded local search.
- No AF3 database directory was found in the bounded local search.
- No real AF3 input JSON jobs or output directories were found in this checkout.

## Result

No real AF3 job was launched. No Phase C audit/export was run on real output.

This is a runtime/acquisition blocker, not a Phase C or Phase D tooling failure.

## Required Unblockers

- AF3 model parameter directory.
- AF3 database root, or valid precomputed data-pipeline inputs.
- Python >=3.12 environment with `lzma` and AF3 dependencies installed.
- One small monomer AF3 input JSON.

After those are present, run one small AF3 job and then run:

```bash
PYTHONPATH=src python tools/audit_af3_output_geometry_readiness.py \
  --af3-output-root /data/af3_single_run/<system_id> \
  --output-root artifacts/nmr_geometry_readiness/phase_f_single_run
```

```bash
PYTHONPATH=src python tools/export_af3_geometry_bundle.py \
  --af3-output-dir /data/af3_single_run/<system_id> \
  --output-dir artifacts/nmr_geometry_readiness/phase_f_single_run_bundle \
  --system-id <system_id> \
  --emit-plots
```
