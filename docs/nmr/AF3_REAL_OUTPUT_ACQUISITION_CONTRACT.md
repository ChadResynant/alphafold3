# AF3 Real Output Acquisition Contract

This contract defines how real AlphaFold3 jobs must be generated so the Phase D
real geometry-readiness survey can produce a scientific cohort result.

Scope is AF3-native output acquisition only. Do not include NMR residuals,
reference structures, reference RMSD, reference mappings, Predyctor, BPHON, or
Valydator artifacts.

## Scientific Target

The survey should answer whether real AF3 ranked outputs are mapping-grade,
geometrically complete, and near-duplicate versus ensemble-like.

The acquisition must preserve enough AF3-native information to test:

- atom identity integrity
- coordinate completeness
- confidence linkage
- ranked-model diversity
- confidence versus real geometry variance

## Required AF3 Job Layout

Each system must be a standard AF3 output directory compatible with
`docs/output.md`:

```text
<af3_output_root>/<system_id>/
  <system_id>_data.json
  <system_id>_ranking_scores.csv
  <system_id>_model.cif
  <system_id>_confidences.json
  <system_id>_summary_confidences.json
  seed-<seed>_sample-<sample>/
    <system_id>_seed-<seed>_sample-<sample>_model.cif
    <system_id>_seed-<seed>_sample-<sample>_confidences.json
    <system_id>_seed-<seed>_sample-<sample>_summary_confidences.json
```

`.cif.zst` and `_confidences.json.zst` are allowed for large files. Summary
confidence JSON must remain present.

The Phase D survey treats each `<system_id>` directory as one AF3 job. It does
not infer jobs from unrelated filenames.

## Minimum Real Cohort

Acquire enough real AF3 outputs to support cohort-level interpretation:

- at least 10 single-chain protein systems
- at least 5 multichain or complex systems, if available
- at least 3 nucleic-acid/protein or ligand-containing systems, if supported
- at least 3 ranked models per system

The smallest useful unblocker is five real AF3 output directories. Fewer than
five real outputs is a no-go for scientific interpretation.

## Required Metadata

For each system, record an operator manifest row with:

```text
system_id
input_json_path
input_json_sha256
af3_output_dir
ranking_scores_path
model_rank
seed
sample
model_cif_path
confidence_json_path
summary_confidence_json_path
```

Rules:

- `system_id` must be stable and must match the AF3 output directory name.
- `input_json_sha256` is computed from the exact AF3 input JSON submitted.
- `model_rank` must be derived from `<system_id>_ranking_scores.csv`.
- model CIF and confidence JSON paths must point to per-seed/per-sample outputs,
  not only the top-ranked root copy.
- Do not add reference PDB/mmCIF paths or NMR residual paths to this manifest.

## AF3 Run Requirements

Use at least three ranked model candidates per system. With the current AF3
driver, that means using enough seeds and/or diffusion samples to produce at
least three rows in `<system_id>_ranking_scores.csv`.

Example:

```bash
python run_alphafold.py \
  --json_path=/data/af3_inputs/<system_id>.json \
  --output_dir=/data/af3_real_outputs \
  --model_dir=/data/af3_models \
  --db_dir=/data/af3_databases \
  --num_diffusion_samples=5 \
  --num_seeds=1 \
  --run_data_pipeline=true \
  --run_inference=true
```

If using precomputed data-pipeline outputs, keep the same output layout:

```bash
python run_alphafold.py \
  --json_path=/data/af3_inputs/<system_id>.json \
  --output_dir=/data/af3_real_outputs \
  --model_dir=/data/af3_models \
  --db_dir=/data/af3_databases \
  --num_diffusion_samples=5 \
  --num_seeds=1 \
  --run_data_pipeline=false \
  --run_inference=true
```

## Survey Commands

Run the Phase D real-output survey:

```bash
PYTHONPATH=src python tools/survey_real_af3_geometry_readiness.py \
  --search-root /data/af3_real_outputs \
  --output-root artifacts/nmr_geometry_readiness/phase_d_real_survey
```

If previous synthetic or survey artifacts live under the same root, exclude
them:

```bash
PYTHONPATH=src python tools/survey_real_af3_geometry_readiness.py \
  --search-root /data/af3_real_outputs \
  --exclude-path artifacts/nmr_geometry_readiness \
  --output-root artifacts/nmr_geometry_readiness/phase_d_real_survey
```

Optionally export a single system bundle for inspection:

```bash
PYTHONPATH=src python tools/export_af3_geometry_bundle.py \
  --af3-output-dir /data/af3_real_outputs/<system_id> \
  --output-dir artifacts/nmr_geometry_readiness/manual_bundle/<system_id> \
  --system-id <system_id> \
  --emit-plots
```

## Required Survey Outputs

The Phase D survey must produce:

```text
artifacts/nmr_geometry_readiness/phase_d_real_survey/
  REPORT.md
  metrics.json
  metrics_cohort.csv
  blocked_reasons.csv
  diversity_grade_summary.csv
  confidence_geometry_summary.csv
```

If real AF3 output directories are found, it may also produce:

```text
artifacts/nmr_geometry_readiness/phase_d_real_survey/bundles/<system_id>/
```

## No-Go Conditions

Stop and report no-go if any of these are true:

- fewer than 5 real AF3 output directories are available
- any counted system has fewer than 3 ranked models
- per-rank model CIF files are missing
- per-rank confidence JSON files are missing
- `<system_id>_ranking_scores.csv` is missing or unparsable
- atom identity is unstable: duplicate canonical atom identities, missing
  required `_atom_site` columns, nontrivial altloc ambiguity, or invalid
  coordinates

Do not repair by sequence alignment, coordinate matching, residue renumbering,
filename heuristics, confidence metadata, or reference structures.

## Interpretation Gate

Once the survey passes acquisition no-go checks, interpret the cohort by:

- `readiness_verdict_counts`
- `geometry_grade_counts`
- `diversity_grade_counts`
- `confidence_geometry_agreement_counts`
- common blocked reasons

The immediate scientific readout is whether real AF3 ranked outputs are mostly:

- `near_duplicate_models`: AF3 behaves as a high-quality single-state predictor
- `modest_local_diversity`: AF3 carries limited local uncertainty
- `ensemble_like`: AF3 ranked outputs approximate a structural ensemble
