# AF3 Geometry Readiness Analysis Plan

## Core Scientific Question

Can AF3 be treated as a governed, atom-resolved, ensemble-like geometry source
for downstream physical modeling?

The analysis is AF3-native. It does not evaluate NMR residuals, compare against
reference structures, compute reference RMSD, or infer chain/residue mappings.

## Figure Set

1. AF3 output completeness/readiness
   - extraction flow
   - readiness verdict counts
   - blocked reason counts
2. Atom/residue/chain identity stability
   - chain, residue, atom counts
   - duplicate canonical atom keys
   - auth-vs-label discrepancies
   - residue numbering policy
3. Ranked-model diversity
   - rank-1 vs rank-k backbone RMSD
   - pairwise backbone RMSD distribution
   - per-residue coordinate variance
4. Local diversity by structural class
   - backbone variance
   - aromatic centroid and plane variance
   - methyl carbon variance
   - amide vector variance only with explicit H/HN atoms
5. Confidence versus real geometric diversity
   - local pLDDT versus coordinate variance
   - confidence/geometry agreement grade
6. Sidechain microstate behavior
   - aromatic and methyl local diversity summaries

## Thresholds

- `mapping_grade` requires no duplicate atom keys, no invalid coordinates,
  confidence atom-count match, and all required `_atom_site` columns.
- Diversity grading requires at least two ranked models, 20 shared backbone
  atoms, and five shared residues.
- `ensemble_like` requires median pairwise backbone RMSD at least `0.5 A` or p95
  local residue variance at least `0.25 A^2`.
- `near_duplicate_models` requires median pairwise backbone RMSD below `0.1 A`
  and p95 local residue variance below `0.05 A^2`.
- Confidence/geometry agreement requires at least 30 residues with both pLDDT
  and coordinate variance.

## Interpretation

Positive outcome: AF3 ranked outputs contain meaningful local structural
diversity and can be packaged as deterministic geometry artifacts.

Negative outcome: AF3 outputs are high-quality single-state predictions or are
not mapping-grade without additional curation.

Both outcomes are useful because they separate coordinate usability, identity
stability, confidence linkage, and ensemble-like diversity.
