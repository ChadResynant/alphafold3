# AF3 Runtime Unblock Checklist

Goal: unblock **Phase F** (one real AF3 job) and therefore unblock **Phase D** (real
geometry-readiness survey) by ensuring the AF3 runtime environment can produce a
single output directory that our Phase C audit/export can consume.

Scope: alphafold3 repo only. No downstream NMR evaluation. No reference mapping.

## 1) Python Environment (Required)

Minimum: **Python >= 3.12** with working `lzma` and AF3 dependencies installed.

Must import:
- `lzma`
- `absl`
- `jax`
- `rdkit`
- `zstandard`

Validation:

```bash
python -c "import lzma, absl, jax, rdkit, zstandard; print('ok')"
```

Repo regression (AF3 runtime smoke test):

```bash
python -m pytest run_alphafold_test.py -q
```

If either command fails, do not proceed to job execution.

## 2) AF3 Model Parameters (Required)

Record:
- params path: `<params_dir>`
- readable: `yes/no`
- checksum/hash (if available): `<hash>`

Note: do not copy/rename files in a way that breaks the expected AF3 directory
layout for parameters.

## 3) AF3 Databases / Pipeline Inputs (Required)

Record:
- db path: `<db_dir>`
- readable: `yes/no`
- minimum DB subset sufficient for a small monomer: `<notes>`

If using precomputed pipeline inputs instead of full databases, record:
- pipeline inputs path: `<pipeline_inputs_dir>`
- provenance / job id: `<id>`

## 4) One Small Monomer Input JSON (Required)

Use the simplest “golden path” first:
- single-chain monomer
- ~100-200 residues
- no ligands, no nucleic acids

Record:
- `system_id`: `<system_id>`
- sequence length: `<n_residues>`
- input JSON path: `<input.json>`
- intended output dir: `<output_dir>`

## 5) First Real Run (Target)

Run:

```bash
python run_alphafold.py \
  --json_path <input.json> \
  --model_dir <params_dir> \
  --db_dir <db_dir> \
  --output_dir <output_dir>
```

Minimum expected output layout (per system):
- at least 3 ranked model CIFs (e.g. `rank_001.cif`, `rank_002.cif`, `rank_003.cif`)
- confidence JSON per ranked model, plus any AF3 summary confidence artifact if
  emitted by the runtime
- ranking metadata if emitted by the runtime

## 6) Phase C Audit/Export on the Real Output (Required)

Audit:

```bash
PYTHONPATH=src python tools/audit_af3_output_geometry_readiness.py \
  --af3-output-root <output_dir>/<system_id> \
  --output-root artifacts/nmr_geometry_readiness/phase_f_single_run_real \
  --emit-plots
```

Export bundle:

```bash
PYTHONPATH=src python tools/export_af3_geometry_bundle.py \
  --af3-output-dir <output_dir>/<system_id> \
  --output-dir artifacts/nmr_geometry_readiness/phase_f_single_run_real_bundle \
  --system-id <system_id> \
  --emit-plots
```

Success criteria (Phase F complete):
- `readiness_verdict` is not blocked
- `identity_grade` is `stable` (mapping-grade)
- `geometry_grade` is `complete`
- `diversity_grade` is one of:
  - `near_duplicate_models`
  - `modest_local_diversity`
  - `ensemble_like`
  - `not_applicable` (only if gating prevents grading)
- `confidence_geometry_agreement` is computed or `not_applicable` (per gating)

Do not infer or “repair” mappings, identity, or geometry. Fail closed.

## 7) Phase D Real Survey (After One Successful Run)

Once at least one real AF3 output exists:

```bash
PYTHONPATH=src python tools/survey_real_af3_geometry_readiness.py \
  --search-root <output_dir> \
  --output-root artifacts/nmr_geometry_readiness/phase_g_first_real_output
```

Stop after producing the report surfaces. Do not add downstream evaluation here.

