# Next Phase: First Real AF3 Output (Stop After One)

This repo already contains:
- Phase C: geometry-readiness audit + bundle export + diversity/confidence metrics
- Phase D: real-output survey (blocked until real outputs exist)
- Phase E: real-output acquisition contract

The next work is **execution**, not more instrumentation.

## Goal

Produce exactly **one** real AF3 output directory that:
1. is consumable by Phase C tools
2. yields a real-data readiness classification

Then immediately run the Phase D survey over the output root and stop.

No cohort expansion and no downstream evaluation in this repo.

## Prerequisites (Hard Gate)

Follow the one-page checklist:
- [AF3_RUNTIME_UNBLOCK_CHECKLIST.md](AF3_RUNTIME_UNBLOCK_CHECKLIST.md)

Do not proceed unless:
- `python -c "import lzma, absl, jax, rdkit, zstandard; print('ok')"` succeeds
- `python -m pytest run_alphafold_test.py -q` succeeds
- AF3 model parameters directory is present
- AF3 database root (or valid precomputed pipeline inputs) is present

## Execution Steps

1. Pick a small monomer (~100-200 aa), no ligands, no nucleic acids.
2. Create one AF3 input JSON with a stable `system_id`.
3. Run AF3 to generate >=3 ranked models:

```bash
python run_alphafold.py \
  --json_path <input.json> \
  --model_dir <params_dir> \
  --db_dir <db_dir> \
  --output_dir <output_dir>
```

4. Run Phase C audit on the **real** output directory:

```bash
PYTHONPATH=src python tools/audit_af3_output_geometry_readiness.py \
  --af3-output-root <output_dir>/<system_id> \
  --output-root artifacts/nmr_geometry_readiness/phase_f_single_run_real \
  --emit-plots
```

5. Run Phase C bundle export on the **real** output directory:

```bash
PYTHONPATH=src python tools/export_af3_geometry_bundle.py \
  --af3-output-dir <output_dir>/<system_id> \
  --output-dir artifacts/nmr_geometry_readiness/phase_f_single_run_real_bundle \
  --system-id <system_id> \
  --emit-plots
```

6. Run Phase D survey over the output root (now unblocked by at least one job):

```bash
PYTHONPATH=src python tools/survey_real_af3_geometry_readiness.py \
  --search-root <output_dir> \
  --output-root artifacts/nmr_geometry_readiness/phase_g_first_real_output
```

## Stop Conditions

Stop the effort for the day after either:
- one successful real job + Phase C audit/export + Phase D survey report surfaces exist, or
- a definitive runtime blocker is identified (environment, params, db, input format).

Do not start a 5-10 job cohort until the “one real output” path is clean.

## Interpreting the First Real Result

Primary outputs to inspect:
- `REPORT.md`
- `metrics.json`
- `geometry_diversity_summary.csv`
- `confidence_geometry_correlation.csv`
- `figures/`

Decision questions:
1. mapping-grade identity?
2. geometry complete/usable?
3. ensemble-like vs near-duplicate behavior?
4. confidence aligned with real coordinate variance (or not applicable by gating)?

If the first real run is `not_mapping_grade` or `unusable`, that is a legitimate
AF3-output characterization result and should be reported as such (fail-closed,
no repair/inference).

