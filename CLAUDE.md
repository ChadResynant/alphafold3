# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

## Architecture Overview

AlphaFold3 is a multi-stage structure prediction system. The pipeline is split into:

**1. Data Pipeline** (`alphafold3/data/`): Featurization and preprocessing
   - `pipeline.py` — orchestrates genetic search (MSA, templates), controls `--run_data_pipeline` flag
   - `featurisation.py` — converts input sequences/structures into numeric features
   - `msa.py` — Multiple Sequence Alignment via HMMer suite (jackhmmer, hmmsearch, hmmalign)
   - `templates.py` — Structural template search and alignment
   - C++ extensions (`data/cpp/`) wrap PDB CIF parsing (libcifpp), DSSP (secondary structure)

**2. Model Inference** (`alphafold3/model/`): JAX-based neural network
   - `model.py` — main entry point; loads weights, runs forward pass
   - `network/` — transformer blocks, confidence head, structure decoder
   - `components/` — attention, embedding, MLPs; JAX utilities
   - `atom_layout/` — maps atoms to canonical positions/ordering
   - `features.py` / `feat_batch.py` — feature representation and batching

**3. Structure Post-Processing** (`alphafold3/structure/`): Geometry refinement
   - C++ component wrapping libcifpp, Einstein summation optimization
   - Outputs relaxed structures, confidence metrics

**4. Constants & Utilities**
   - `constants/` — chemical properties, atom vocabularies
   - `common/` — I/O (folding_input JSON parsing), logging, resource management
   - `parsers/` — C++ wrappers for structure/metadata parsing (pybind11 bindings)

**Key insight:** The JAX model (`model.py`) is stateless and operates on numpy arrays. Weights are loaded separately via `model/params.py`. C++ bindings (pybind11) handle I/O-heavy tasks like PDB parsing and secondary structure assignment.

## Development Commands

**Python version:** AF3 requires Python ≥3.12 (see `pyproject.toml`)

```bash
# Pre-flight checks (always first)
bash ~/repos/build-lib/preflight/surface_all.sh .

# Install in editable mode (builds C++ extensions via scikit-build-core)
pip install -e .

# Install with test dependencies
pip install -e ".[dev]"

# Build without installing (for development)
pip install -e . --no-build-isolation --verbose

# Build intermediate data (chemical component pickles for constants)
build_data

# Run integration tests (end-to-end pipeline)
# Uses absltest framework (from absl-py); run directly with Python
python run_alphafold_test.py

# Run data pipeline tests only
python run_alphafold_data_test.py

# Run tests with verbose output
python run_alphafold_test.py --alsologtostderr --v=2

# Run a specific test case
python run_alphafold_test.py AlphaFoldTest.test_case_name
```

## Testing Framework

AF3 uses **absltest** (Google's testing framework from `absl-py`), not pytest. The main test suites are:

- `run_alphafold_test.py` — Full integration test (end-to-end pipeline with model inference). Requires GPU, model weights, and all databases. Runs the complete structure prediction workflow.
- `run_alphafold_data_test.py` — Data pipeline tests (featurization, MSA, template search). Can run on CPU but requires HMMER binaries (jackhmmer, hmmsearch, etc.).

Both test files use `absltest` and are run directly with Python, not via pytest. Tests use temporary directories managed by the test framework (accessible via `absltest.TEST_TMPDIR`).

**Note:** These are large integration tests. For rapid iteration on specific modules, write unit tests inline using Python's `unittest` framework.

## Running Structure Predictions

```bash
# Basic prediction (requires model weights + databases)
python run_alphafold.py \
  --json_path=fold_input.json \
  --model_dir=$RESYNANT_DATA/alphafold3/models \
  --output_dir=$RESYNANT_ARTIFACTS/alphafold3/predictions

# Skip data pipeline (use pre-computed features)
python run_alphafold.py \
  --json_path=fold_input.json \
  --model_dir=$RESYNANT_DATA/alphafold3/models \
  --output_dir=$RESYNANT_ARTIFACTS/alphafold3/predictions \
  --run_data_pipeline=false

# Skip inference (data pipeline only)
python run_alphafold.py \
  --json_path=fold_input.json \
  --output_dir=$RESYNANT_ARTIFACTS/alphafold3/features \
  --run_inference=false

# See all flags
python run_alphafold.py --help
```

## Key Files & Entry Points

| File | Purpose |
|------|---------|
| `run_alphafold.py` | Main CLI script; parses JSON input, orchestrates pipeline |
| `src/alphafold3/model/model.py` | JAX inference; loads weights, runs network forward pass |
| `src/alphafold3/data/pipeline.py` | Genetic search orchestration (MSA, templates) |
| `src/alphafold3/data/featurisation.py` | Input → features conversion |
| `src/alphafold3/common/folding_input.py` | JSON input parsing |
| `src/alphafold3/model/params.py` | Model weight loading |
| `CMakeLists.txt` | C++ build config (libcifpp, DSSP, pybind11 bindings) |
| `pyproject.toml` | Python build config (scikit-build-core) |

## Build System (Scikit-Build + CMake)

AF3 uses **scikit-build-core** to integrate CMake with Python packaging. The build process:

1. CMake fetches C++ dependencies (abseil, pybind11, libcifpp, DSSP, etc.)
2. Compiles C++ extensions (`alphafold3.cpp` bindings for structure I/O, secondary structure)
3. Python packaging wraps the compiled `.so` files
4. Post-install: `build_data` generates chemical component pickles from libcifpp data

**Build requirements:**
- Python ≥3.12
- cmake ≥3.28, ninja, C++20 compiler
- macOS: `xcode-select --install`
- Linux: `apt-get install cmake ninja-build`
- HMMER suite (jackhmmer, hmmsearch, hmmalign, hmmbuild) for data pipeline tests

**Troubleshooting:**
- If build fails, check `CMakeLists.txt` line 15 — FetchContent Git config workaround (fixed in recent CMake)
- `build_data` requires `LIBCIFPP_DATA_DIR` or searches standard install paths for `libcifpp/components.cif`
- Force rebuild: `pip install -e . --no-cache-dir --force-reinstall`

## RESYNANT Integration

**Model Weights & Databases:**
- Store model parameters in `$RESYNANT_DATA/alphafold3/models/`
- Store genetic databases (BFD, UniRef, etc.) in `$RESYNANT_DATA/alphafold3/databases/`
- Download via `fetch_databases.sh` — modify to save to RAID instead of local disk

**Structure Predictions (Artifacts):**
- Write predictions to `$RESYNANT_ARTIFACTS/alphafold3/{experiment_name}/`
- Include provenance: `source ~/repos/scripts/run-provenance.sh alphafold3 <run_id>`
- Immutable after promotion (add SHA manifest)

**No Hardcoded Paths:**
- Use `$RESYNANT_DATA`, `$RESYNANT_ARTIFACTS` environment variables
- Reference: [`governance/contracts/PATH_RESOLUTION_CONTRACT.md`](../governance/contracts/PATH_RESOLUTION_CONTRACT.md)

## Cross-Repo Integration

AF3 can feed predictions to shift prediction pipelines (shiftx4proteins, quantumchemistry, predyctor):
- Predictions are returned as mmCIF files (with confidence metadata)
- Atom identities follow [`governance/contracts/ATOM_IDENTITY_CONTRACT.md`](../governance/contracts/ATOM_IDENTITY_CONTRACT.md)
- See `alphafold3/structure/` for structure representation

## Key Policies & Contracts

- **Path Resolution:** [`governance/contracts/PATH_RESOLUTION_CONTRACT.md`](../governance/contracts/PATH_RESOLUTION_CONTRACT.md)
- **Atom Identity:** [`governance/contracts/ATOM_IDENTITY_CONTRACT.md`](../governance/contracts/ATOM_IDENTITY_CONTRACT.md)
- **Build Protocols:** [`build-lib/docs/BUILD_PROTOCOLS.md`](../build-lib/docs/BUILD_PROTOCOLS.md)

## Development Patterns

**Rapid iteration:** For quick testing of individual modules (e.g., featurization, scoring), use `--run_inference=false` or `--run_data_pipeline=false` to skip expensive stages.

**Debugging data pipeline:** Enable logging with `--alsologtostderr --v=2` when running tests or predictions to see detailed messages from the data pipeline stages (MSA, templates, featurization).

**JAX/GPU development:**
- JAX installation is finicky; follow [`docs/installation.md`](docs/installation.md)
- Verify GPU: `python -c "import jax; print(jax.devices())"`
- Force CPU-only: `JAX_PLATFORMS=cpu`
- Profile GPU memory: set `JAX_CHECK_TRACER_LEAKS=1` to catch traced values leaking

**Memory management:**
- Genetic search (data pipeline) requires ~64 GB RAM for long protein sequences
- Use `--run_data_pipeline=false` for quick inference-only testing
- Features cache during `--run_inference=false` mode; reuse cached features to avoid re-searching

**C++ compilation issues:**
- Check actual compile commands: `cmake --build . -v`
- FetchContent Git issues on newer Git: Fixed by `GIT_CONFIG_PARAMETERS` in CMakeLists.txt line 15
- libcifpp/DSSP may fail on some systems; check `LIBCIFPP_DATA_DIR` environment variable for `build_data`

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for submission guidelines. Key points:
- Small patches (bug fixes, documentation) are welcome
- AI-generated code must be transparently labeled, reviewed, and tested
- All contributions require Contributor License Agreement (CLA)
- PRs require manual review

## Resources

- [Installation docs](docs/installation.md)
- [Input format spec](docs/input.md)
- [Output format spec](docs/output.md)
- [Known issues](docs/known_issues.md)
- [Performance guide](docs/performance.md)
- Upstream: [google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3)
