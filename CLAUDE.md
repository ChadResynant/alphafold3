# CLAUDE.md — alphafold

**AlphaFold3** — Protein structure prediction from Google DeepMind.

This is a clone of [`google-deepmind/alphafold3`](https://github.com/google-deepmind/alphafold3) (upstream: origin), mirrored to Falcon RAID (lan remote).

## Governance

**Before any work:** read [`~/repos/governance/INDEX.md`](../governance/INDEX.md)

All cross-repo contracts, policies, and CI gates in `~/repos/governance/` apply here.

## Remotes

- **origin** → `https://github.com/google-deepmind/alphafold3.git` (upstream Google)
- **lan** → `/Volumes/falcon-raid/resynant/git/alphafold3.git` (RAID mirror, macOS)

To sync with upstream:
```bash
git fetch origin
git merge origin/main  # or rebase if tracking upstream commits
```

## Quick Reference

```bash
# Pre-flight checks
bash ~/repos/build-lib/preflight/surface_all.sh .

# Install AF3 and dependencies
pip install -e .
# See: docs/install.md

# Run structure prediction
python run_alphafold.py --input_dir=<pdb|seq> --output_dir=<dir>

# Run tests
python run_alphafold_test.py
```

## Key Policies

- **Path Resolution:** [`governance/contracts/PATH_RESOLUTION_CONTRACT.md`](../governance/contracts/PATH_RESOLUTION_CONTRACT.md) 
  - All data must use `$RESYNANT_DATA`, `$RESYNANT_ARTIFACTS`, etc. — never hardcoded paths
- **Atom Identity:** [`governance/contracts/ATOM_IDENTITY_CONTRACT.md`](../governance/contracts/ATOM_IDENTITY_CONTRACT.md)
  - Canonical join key: `(system_id, chain_id, residue_number, atom_name)`. Primes are identity.
- **Build Protocols:** [`build-lib/docs/BUILD_PROTOCOLS.md`](../build-lib/docs/BUILD_PROTOCOLS.md)

## Integration Notes

- **Model Weights:** Downloaded by `fetch_databases.sh`. Store on RAID under `$RESYNANT_DATA/alphafold3/`.
- **Artifacts:** Structure predictions and feature caches → `$RESYNANT_ARTIFACTS/alphafold3/`
- **No Hardcoded Paths:** AF3 integration must not assume `/home/user/` or machine-specific paths.

## File Structure

```
alphafold3/
├── run_alphafold.py           # Main entry point
├── fetch_databases.sh         # Download model weights + databases
├── docs/                      # Documentation
├── docker/                    # Docker build configs
├── legal/                     # License and usage terms
└── alphafold3/                # Python package
    ├── models/                # Model definitions
    ├── data/                  # Data processing
    └── ...
```

## Next Steps

1. Run `fetch_databases.sh` to download model weights
2. Test with `python run_alphafold_test.py`
3. Integrate into NMR shift prediction pipelines (shiftx4proteins, quantumchemistry, etc.)
4. Configure data paths to use RAID mirrors (`$RESYNANT_DATA`, `$RESYNANT_ARTIFACTS`)
