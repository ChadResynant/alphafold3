# AF3 Real Geometry Readiness Survey

Phase D runs the Phase C AF3-native readiness audit over real AF3 output
directories and summarizes the cohort. It does not read NMR residuals, does not
use reference structures, does not compute reference RMSD, and does not infer or
repair mappings.

## Command

```bash
python tools/survey_real_af3_geometry_readiness.py \
  --search-root <real_af3_outputs_root> \
  --output-root artifacts/nmr_geometry_readiness/phase_d_real_survey
```

Use `--exclude-path` for generated synthetic artifacts or previous survey output
roots that should not be treated as real AF3 output.

## Outputs

The survey writes:

- `REPORT.md`
- `metrics.json`
- `metrics_cohort.csv`
- `blocked_reasons.csv`
- `diversity_grade_summary.csv`
- `confidence_geometry_summary.csv`
- `bundles/<system_id>/` for each valid discovered AF3 output directory

If no real AF3 output directories are found, the survey emits
`BLOCKED_NO_REAL_AF3_OUTPUTS` and stops without inferring substitutes.

## Cohort Questions

The report answers:

- how many outputs are `mapping_grade`, `mapping_repairable`, or
  `not_mapping_grade`
- how many are `complete`, `partially_complete`, or `unusable`
- how many ranked-model sets are `near_duplicate_models`,
  `modest_local_diversity`, `ensemble_like`, or `not_applicable`
- whether confidence is `aligned`, `weak`, `decoupled`, or `not_applicable`
- which blocked reasons are most common
