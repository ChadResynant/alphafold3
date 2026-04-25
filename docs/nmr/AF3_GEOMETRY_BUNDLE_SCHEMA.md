# AF3 Geometry Bundle Schema

The AF3 geometry bundle is an AlphaFold3-owned producer artifact for deterministic
atom-level geometry export. It does not contain NMR residuals, reference
structures, reference RMSD, or inferred mappings.

## Layout

```text
af3_geometry_bundle/
  bundle_manifest.json
  af3_models/
    rank_001.cif
    rank_002.cif
  confidence/
    rank_001_confidence.json
    rank_001_summary_confidence.json
  atom_identity_table.parquet
  residue_identity_table.parquet
  chain_identity_table.csv
  model_ranking_table.csv
  geometry_diversity_table.csv
  local_class_diversity.csv
  pairwise_model_rmsd.csv
  per_residue_variance.csv
  confidence_geometry_correlation.csv
  REPORT.md
```

## Manifest Guardrails

`bundle_manifest.json` must include:

```json
{
  "artifact_class": "af3_geometry_bundle",
  "schema_version": "af3.geometry_bundle.v1",
  "producer_repo": "alphafold3",
  "promotion_allowed": false,
  "training_surface_allowed": false,
  "contains_nmr_residuals": false,
  "contains_reference_structure": false,
  "identity_policy": "af3_native_explicit_only",
  "mapping_inference_allowed": false,
  "altloc_policy": "single_conformer_only"
}
```

## Canonical Atom Identity

Canonical identity is:

```text
system_id, model_rank, chain_id, residue_number, insertion_code, atom_name, element
```

where:

- `chain_id` comes from `_atom_site.auth_asym_id`
- `residue_number` comes from `_atom_site.auth_seq_id`
- `atom_name` comes from `_atom_site.label_atom_id`
- `element` comes from `_atom_site.type_symbol`
- insertion code values `.`, `?`, and missing are normalized to the empty string

The bundle retains label IDs separately for auditability.

## Fail-Closed Rules

The exporter returns `not_mapping_grade` for:

- missing required `_atom_site` columns
- duplicate canonical atom identities
- invalid or non-finite coordinates
- confidence atom-count mismatch
- disallowed or ambiguous altlocs

Altloc policy allows only empty, `.`, `?`, or a single `A` conformer.
