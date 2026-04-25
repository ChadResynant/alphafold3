# AF3 Output Geometry Readiness

The readiness audit answers whether AF3 outputs can be consumed as deterministic
atom-resolved geometry artifacts. It is an AF3-native producer check, not a
downstream NMR or reference-structure evaluation.

## Command

```bash
python tools/audit_af3_output_geometry_readiness.py \
  --af3-output-root <output_root> \
  --output-root artifacts/nmr_geometry_readiness/<run_id> \
  --emit-plots
```

## Verdicts

- `mapping_grade`: atom identity, coordinates, and confidence linkage pass the
  fail-closed checks.
- `mapping_repairable`: only repairable packaging blockers are present.
- `not_mapping_grade`: hard identity, coordinate, altloc, or confidence-linkage
  blockers are present.

## Grades

- `identity_grade = stable | repairable | unstable`
- `geometry_grade = complete | partially_complete | unusable`
- `diversity_grade = near_duplicate_models | modest_local_diversity | ensemble_like | not_applicable`
- `confidence_geometry_agreement = aligned | weak | decoupled | not_applicable`

## Output Tables

The audit writes deterministic CSV/JSON outputs:

- `bundle_readiness_summary.csv`
- `identity_metrics.csv`
- `geometry_completeness.csv`
- `geometry_diversity_summary.csv`
- `pairwise_model_rmsd.csv`
- `per_residue_variance.csv`
- `confidence_geometry_correlation.csv`
- `local_class_diversity.csv`
- `readiness.json`
- `REPORT.md`
- `figures/` when `--emit-plots` is enabled

No reference mapping is inferred. No NMR residuals are read or written.
