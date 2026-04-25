# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0.

"""AF3-native geometry bundle schema and mmCIF identity extraction.

This module is intentionally self-contained and does not import the compiled
AlphaFold mmCIF parser. The geometry-readiness tools must be able to inspect
already-written AF3 outputs in source checkouts where the C++ extension has not
been built.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
import csv
import dataclasses
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Any

import zstandard


SCHEMA_VERSION = "af3.geometry_bundle.v1"
ARTIFACT_CLASS = "af3_geometry_bundle"
IDENTITY_POLICY = "af3_native_explicit_only"
ALTLOC_POLICY = "single_conformer_only"

READINESS_MAPPING_GRADE = "mapping_grade"
READINESS_MAPPING_REPAIRABLE = "mapping_repairable"
READINESS_NOT_MAPPING_GRADE = "not_mapping_grade"

IDENTITY_STABLE = "stable"
IDENTITY_REPAIRABLE = "repairable"
IDENTITY_UNSTABLE = "unstable"

GEOMETRY_COMPLETE = "complete"
GEOMETRY_PARTIALLY_COMPLETE = "partially_complete"
GEOMETRY_UNUSABLE = "unusable"

RESIDUE_NUMBERING_AUTH_PRESERVED = "auth_seq_id_preserved"
RESIDUE_NUMBERING_CONTAINS_GAPS = "contains_gaps"
RESIDUE_NUMBERING_NON_MONOTONIC = "non_monotonic"

REQUIRED_ATOM_SITE_COLUMNS = (
    "_atom_site.auth_asym_id",
    "_atom_site.auth_seq_id",
    "_atom_site.label_atom_id",
    "_atom_site.type_symbol",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
)

OPTIONAL_ATOM_SITE_COLUMNS = (
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_alt_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.pdbx_PDB_model_num",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.group_PDB",
)

ALLOWED_ALTLOC_VALUES = {"", ".", "?", "A"}
CANONICAL_IDENTITY_FIELDS = (
    "system_id",
    "model_rank",
    "chain_id",
    "residue_number",
    "insertion_code",
    "atom_name",
    "element",
)

BLOCKED_MISSING_ATOM_SITE_COLUMNS = "missing_atom_site_columns"
BLOCKED_DUPLICATE_ATOM_IDENTITIES = "duplicate_atom_identities"
BLOCKED_CONFIDENCE_ATOM_COUNT_MISMATCH = "confidence_atom_count_mismatch"
BLOCKED_NONTRIVIAL_ALTLOC = "nontrivial_altloc"
BLOCKED_INVALID_COORDINATES = "invalid_coordinates"
BLOCKED_MISSING_RANKED_MODELS = "missing_ranked_models"

_REPAIRABLE_BLOCKERS = frozenset({BLOCKED_MISSING_RANKED_MODELS})


@dataclasses.dataclass(frozen=True, slots=True)
class SourceModel:
  """Paths and ranking metadata for a single AF3 model output."""

  rank_index: int
  model_rank: str
  seed: str
  sample: str
  ranking_score: float | None
  model_path: Path
  confidence_path: Path | None
  summary_confidence_path: Path | None


@dataclasses.dataclass(frozen=True, slots=True)
class ModelExtraction:
  """Extracted identity rows and validation status for one model."""

  source_model: SourceModel
  atom_rows: list[dict[str, Any]]
  residue_rows: list[dict[str, Any]]
  chain_rows: list[dict[str, Any]]
  blocked_reasons: tuple[str, ...]
  metrics: dict[str, Any]


@dataclasses.dataclass(frozen=True, slots=True)
class BundleTables:
  """All deterministic tables emitted by an AF3 geometry bundle."""

  atom_rows: list[dict[str, Any]]
  residue_rows: list[dict[str, Any]]
  chain_rows: list[dict[str, Any]]
  ranking_rows: list[dict[str, Any]]
  extraction_metrics: list[dict[str, Any]]
  blocked_reasons: tuple[str, ...]
  atom_table_hash: str
  residue_numbering_policy: str
  readiness_verdict: str
  identity_grade: str
  geometry_grade: str


def read_text_maybe_zst(path: Path) -> str:
  """Reads UTF-8 text from a plain or zstd-compressed path."""
  if path.suffix == ".zst":
    with zstandard.open(path, "rt", encoding="utf-8") as f:
      return f.read()
  return path.read_text(encoding="utf-8")


def read_json_maybe_zst(path: Path) -> Mapping[str, Any]:
  return json.loads(read_text_maybe_zst(path))


def copy_text_maybe_zst(source: Path, destination: Path) -> None:
  """Copies a text artifact, decompressing zstd sources to stable plain text."""
  destination.parent.mkdir(parents=True, exist_ok=True)
  destination.write_text(read_text_maybe_zst(source), encoding="utf-8")


def stable_json_dump(data: Mapping[str, Any], path: Path) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text(
      json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
  )


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  fieldnames: list[str] = []
  for row in rows:
    for key in row:
      if key not in fieldnames:
        fieldnames.append(key)
  with path.open("w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
      writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_parquet(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
  """Writes Parquet with a clear optional-dependency error."""
  try:
    import pyarrow as pa
    import pyarrow.parquet as pq
  except ModuleNotFoundError as e:
    raise RuntimeError(
        "Parquet bundle export requires the optional AF3 NMR dependency "
        "`pyarrow`. Install the repository with the `nmr` extra."
    ) from e

  path.parent.mkdir(parents=True, exist_ok=True)
  table = pa.Table.from_pylist(list(rows))
  pq.write_table(table, path)


def parse_atom_site_table(mmcif_text: str) -> dict[str, list[str]]:
  """Extracts the `_atom_site` loop from mmCIF text.

  The parser is deliberately narrow: it reads looped `_atom_site` columns and
  tokenizes CIF quoting well enough for AF3 output atom rows. It does not infer
  missing identity fields.
  """
  lines = mmcif_text.splitlines()
  i = 0
  while i < len(lines):
    if lines[i].strip() != "loop_":
      i += 1
      continue
    i += 1
    fields: list[str] = []
    while i < len(lines):
      stripped = lines[i].strip()
      if not stripped:
        i += 1
        continue
      if stripped.startswith("_"):
        fields.append(stripped.split()[0])
        i += 1
        continue
      break
    if not fields or not any(field.startswith("_atom_site.") for field in fields):
      continue
    if not all(field.startswith("_atom_site.") for field in fields):
      continue
    tokens: list[str] = []
    while i < len(lines):
      stripped = lines[i].strip()
      if not stripped:
        i += 1
        continue
      if stripped == "#":
        break
      if stripped == "loop_" or stripped.startswith("data_") or stripped.startswith("_"):
        break
      tokens.extend(_tokenize_cif_line(lines[i]))
      i += 1
    if len(tokens) % len(fields) != 0:
      raise ValueError(
          "atom_site loop token count is not divisible by column count"
      )
    table = {field: [] for field in fields}
    for row_start in range(0, len(tokens), len(fields)):
      row = tokens[row_start : row_start + len(fields)]
      for field, value in zip(fields, row, strict=True):
        table[field].append(value)
    return table
  raise ValueError("mmCIF does not contain a looped _atom_site table")


def _tokenize_cif_line(line: str) -> list[str]:
  tokens: list[str] = []
  current: list[str] = []
  quote: str | None = None
  i = 0
  while i < len(line):
    char = line[i]
    if quote:
      if char == quote and (i + 1 == len(line) or line[i + 1].isspace()):
        tokens.append("".join(current))
        current = []
        quote = None
      else:
        current.append(char)
      i += 1
      continue
    if char.isspace():
      if current:
        tokens.append("".join(current))
        current = []
      i += 1
      continue
    if char in ("'", '"') and not current:
      quote = char
      i += 1
      continue
    current.append(char)
    i += 1
  if quote:
    raise ValueError("unterminated CIF quote in atom_site row")
  if current:
    tokens.append("".join(current))
  return tokens


def extract_model_tables(
    *,
    system_id: str,
    source_model: SourceModel,
) -> ModelExtraction:
  """Extracts deterministic identity rows for one AF3 model."""
  blocked: list[str] = []
  metrics: dict[str, Any] = {
      "system_id": system_id,
      "model_rank": source_model.model_rank,
      "source_model_path": str(source_model.model_path),
  }
  try:
    atom_site = parse_atom_site_table(read_text_maybe_zst(source_model.model_path))
  except Exception as e:  # pylint: disable=broad-exception-caught
    return ModelExtraction(
        source_model=source_model,
        atom_rows=[],
        residue_rows=[],
        chain_rows=[],
        blocked_reasons=(BLOCKED_MISSING_ATOM_SITE_COLUMNS,),
        metrics={
            **metrics,
            "parse_error": str(e),
            "atom_count": 0,
            "required_atom_site_columns_present": False,
        },
    )

  missing_columns = [
      column for column in REQUIRED_ATOM_SITE_COLUMNS if column not in atom_site
  ]
  if missing_columns:
    blocked.append(BLOCKED_MISSING_ATOM_SITE_COLUMNS)
  row_count = _atom_site_row_count(atom_site)
  metrics.update(
      {
          "atom_count": row_count,
          "required_atom_site_columns_present": not missing_columns,
          "missing_atom_site_columns": ";".join(missing_columns),
      }
  )
  if missing_columns:
    return ModelExtraction(
        source_model=source_model,
        atom_rows=[],
        residue_rows=[],
        chain_rows=[],
        blocked_reasons=tuple(sorted(set(blocked))),
        metrics=metrics,
    )

  confidence = _load_confidence(source_model.confidence_path)
  confidence_atom_plddts = confidence.get("atom_plddts") if confidence else None
  confidence_chain_ids = confidence.get("atom_chain_ids") if confidence else None
  confidence_count_matches = True
  if confidence_atom_plddts is not None and len(confidence_atom_plddts) != row_count:
    confidence_count_matches = False
    blocked.append(BLOCKED_CONFIDENCE_ATOM_COUNT_MISMATCH)
  if confidence_chain_ids is not None and len(confidence_chain_ids) != row_count:
    confidence_count_matches = False
    blocked.append(BLOCKED_CONFIDENCE_ATOM_COUNT_MISMATCH)
  metrics["confidence_atom_count_matches"] = confidence_count_matches

  atom_rows: list[dict[str, Any]] = []
  invalid_coordinate_count = 0
  disallowed_altloc_count = 0
  for idx in range(row_count):
    atom_name = _value(atom_site, "_atom_site.label_atom_id", idx)
    element = _value(atom_site, "_atom_site.type_symbol", idx).upper()
    chain_id = _value(atom_site, "_atom_site.auth_asym_id", idx)
    residue_number = _value(atom_site, "_atom_site.auth_seq_id", idx)
    insertion_code = normalize_insertion_code(
        _value(atom_site, "_atom_site.pdbx_PDB_ins_code", idx, default="")
    )
    altloc = normalize_altloc(
        _value(atom_site, "_atom_site.label_alt_id", idx, default="")
    )
    if altloc not in ALLOWED_ALTLOC_VALUES:
      disallowed_altloc_count += 1
    try:
      x = float(_value(atom_site, "_atom_site.Cartn_x", idx))
      y = float(_value(atom_site, "_atom_site.Cartn_y", idx))
      z = float(_value(atom_site, "_atom_site.Cartn_z", idx))
      if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
        raise ValueError("non-finite coordinate")
    except ValueError:
      x = y = z = float("nan")
      invalid_coordinate_count += 1

    row = {
        "system_id": system_id,
        "model_rank": source_model.model_rank,
        "rank_index": source_model.rank_index,
        "atom_index": idx + 1,
        "seed": source_model.seed,
        "sample": source_model.sample,
        "chain_id": chain_id,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
        "atom_name": atom_name,
        "element": element,
        "label_asym_id": _value(atom_site, "_atom_site.label_asym_id", idx, default=""),
        "label_seq_id": _value(atom_site, "_atom_site.label_seq_id", idx, default=""),
        "residue_name": _value(atom_site, "_atom_site.label_comp_id", idx, default=""),
        "model_num": _value(
            atom_site, "_atom_site.pdbx_PDB_model_num", idx, default="1"
        ),
        "altloc": altloc,
        "x": x,
        "y": y,
        "z": z,
        "occupancy": _float_or_none(
            _value(atom_site, "_atom_site.occupancy", idx, default="")
        ),
        "b_factor": _float_or_none(
            _value(atom_site, "_atom_site.B_iso_or_equiv", idx, default="")
        ),
        "group_pdb": _value(atom_site, "_atom_site.group_PDB", idx, default=""),
    }
    if confidence_atom_plddts is not None and idx < len(confidence_atom_plddts):
      row["plddt"] = _float_or_none(confidence_atom_plddts[idx])
    else:
      row["plddt"] = row["b_factor"]
    row["atom_identity_hash"] = atom_identity_hash(row)
    atom_rows.append(row)

  if invalid_coordinate_count:
    blocked.append(BLOCKED_INVALID_COORDINATES)
  if disallowed_altloc_count:
    blocked.append(BLOCKED_NONTRIVIAL_ALTLOC)

  duplicate_count, duplicate_altloc_count = _duplicate_identity_counts(atom_rows)
  if duplicate_count:
    blocked.append(BLOCKED_DUPLICATE_ATOM_IDENTITIES)
  if duplicate_altloc_count:
    blocked.append(BLOCKED_NONTRIVIAL_ALTLOC)

  atom_rows = sorted(atom_rows, key=canonical_atom_sort_key)
  residue_rows = residue_rows_from_atoms(atom_rows)
  chain_rows = chain_rows_from_atoms(atom_rows)
  metrics.update(
      {
          "chain_count": len(chain_rows),
          "residue_count": len(residue_rows),
          "duplicate_atom_identity_count": duplicate_count,
          "duplicate_altloc_identity_count": duplicate_altloc_count,
          "invalid_coordinate_count": invalid_coordinate_count,
          "disallowed_altloc_count": disallowed_altloc_count,
          "auth_vs_label_chain_discrepancies": sum(
              1
              for row in atom_rows
              if row["label_asym_id"] and row["chain_id"] != row["label_asym_id"]
          ),
          "auth_vs_label_residue_discrepancies": sum(
              1
              for row in atom_rows
              if row["label_seq_id"]
              and row["label_seq_id"] != "."
              and row["residue_number"] != row["label_seq_id"]
          ),
      }
  )
  return ModelExtraction(
      source_model=source_model,
      atom_rows=atom_rows,
      residue_rows=residue_rows,
      chain_rows=chain_rows,
      blocked_reasons=tuple(sorted(set(blocked))),
      metrics=metrics,
  )


def build_bundle_tables(
    *,
    system_id: str,
    source_models: Sequence[SourceModel],
) -> BundleTables:
  """Extracts and validates all deterministic bundle tables."""
  if not source_models:
    return BundleTables(
        atom_rows=[],
        residue_rows=[],
        chain_rows=[],
        ranking_rows=[],
        extraction_metrics=[],
        blocked_reasons=(BLOCKED_MISSING_RANKED_MODELS,),
        atom_table_hash=hashlib.sha256(b"").hexdigest(),
        residue_numbering_policy=RESIDUE_NUMBERING_AUTH_PRESERVED,
        readiness_verdict=READINESS_NOT_MAPPING_GRADE,
        identity_grade=IDENTITY_UNSTABLE,
        geometry_grade=GEOMETRY_UNUSABLE,
    )

  extractions = [
      extract_model_tables(system_id=system_id, source_model=source_model)
      for source_model in sorted(source_models, key=lambda m: m.rank_index)
  ]
  atom_rows = sorted(
      [row for extraction in extractions for row in extraction.atom_rows],
      key=canonical_atom_sort_key,
  )
  residue_rows = sorted(
      [row for extraction in extractions for row in extraction.residue_rows],
      key=lambda row: (
          row["system_id"],
          row["model_rank"],
          row["chain_id"],
          _residue_number_sort_value(row["residue_number"]),
          row["insertion_code"],
          row["residue_name"],
      ),
  )
  chain_rows = sorted(
      [row for extraction in extractions for row in extraction.chain_rows],
      key=lambda row: (row["system_id"], row["model_rank"], row["chain_id"]),
  )
  ranking_rows = [
      {
          "system_id": system_id,
          "model_rank": source_model.model_rank,
          "rank_index": source_model.rank_index,
          "seed": source_model.seed,
          "sample": source_model.sample,
          "ranking_score": source_model.ranking_score,
          "model_path": str(source_model.model_path),
          "confidence_path": (
              str(source_model.confidence_path)
              if source_model.confidence_path is not None
              else ""
          ),
          "summary_confidence_path": (
              str(source_model.summary_confidence_path)
              if source_model.summary_confidence_path is not None
              else ""
          ),
      }
      for source_model in sorted(source_models, key=lambda model: model.rank_index)
  ]
  blocked = sorted(
      {
          reason
          for extraction in extractions
          for reason in extraction.blocked_reasons
      }
  )
  policy = residue_numbering_policy(atom_rows)
  atom_hash = atom_table_hash(atom_rows)
  readiness = readiness_verdict(blocked)
  return BundleTables(
      atom_rows=atom_rows,
      residue_rows=residue_rows,
      chain_rows=chain_rows,
      ranking_rows=ranking_rows,
      extraction_metrics=[extraction.metrics for extraction in extractions],
      blocked_reasons=tuple(blocked),
      atom_table_hash=atom_hash,
      residue_numbering_policy=policy,
      readiness_verdict=readiness,
      identity_grade=identity_grade(blocked, policy),
      geometry_grade=geometry_grade(blocked),
  )


def export_bundle_artifacts(
    *,
    bundle_dir: Path,
    system_id: str,
    source_output_dir: Path,
    source_models: Sequence[SourceModel],
    tables: BundleTables,
    diversity_rows: Sequence[Mapping[str, Any]] = (),
    local_class_rows: Sequence[Mapping[str, Any]] = (),
    pairwise_rows: Sequence[Mapping[str, Any]] = (),
    residue_variance_rows: Sequence[Mapping[str, Any]] = (),
    confidence_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
  """Writes the bundle directory and returns the manifest."""
  bundle_dir.mkdir(parents=True, exist_ok=True)
  model_dir = bundle_dir / "af3_models"
  confidence_dir = bundle_dir / "confidence"
  model_dir.mkdir(exist_ok=True)
  confidence_dir.mkdir(exist_ok=True)

  for source_model in sorted(source_models, key=lambda m: m.rank_index):
    copy_text_maybe_zst(
        source_model.model_path,
        model_dir / f"{source_model.model_rank}.cif",
    )
    if source_model.confidence_path is not None:
      copy_text_maybe_zst(
          source_model.confidence_path,
          confidence_dir / f"{source_model.model_rank}_confidence.json",
      )
    if source_model.summary_confidence_path is not None:
      copy_text_maybe_zst(
          source_model.summary_confidence_path,
          confidence_dir / f"{source_model.model_rank}_summary_confidence.json",
      )

  write_parquet(bundle_dir / "atom_identity_table.parquet", tables.atom_rows)
  write_parquet(bundle_dir / "residue_identity_table.parquet", tables.residue_rows)
  write_csv(bundle_dir / "chain_identity_table.csv", tables.chain_rows)
  write_csv(bundle_dir / "model_ranking_table.csv", tables.ranking_rows)
  write_csv(bundle_dir / "geometry_diversity_table.csv", diversity_rows)
  write_csv(bundle_dir / "local_class_diversity.csv", local_class_rows)
  write_csv(bundle_dir / "pairwise_model_rmsd.csv", pairwise_rows)
  write_csv(bundle_dir / "per_residue_variance.csv", residue_variance_rows)
  write_csv(bundle_dir / "confidence_geometry_correlation.csv", confidence_rows)

  manifest = bundle_manifest(
      system_id=system_id,
      source_output_dir=source_output_dir,
      source_models=source_models,
      tables=tables,
  )
  stable_json_dump(manifest, bundle_dir / "bundle_manifest.json")
  return manifest


def bundle_manifest(
    *,
    system_id: str,
    source_output_dir: Path,
    source_models: Sequence[SourceModel],
    tables: BundleTables,
) -> dict[str, Any]:
  return {
      "altloc_policy": ALTLOC_POLICY,
      "artifact_class": ARTIFACT_CLASS,
      "atom_table_hash": tables.atom_table_hash,
      "blocked_reasons": list(tables.blocked_reasons),
      "contains_nmr_residuals": False,
      "contains_reference_structure": False,
      "geometry_grade": tables.geometry_grade,
      "identity_grade": tables.identity_grade,
      "identity_policy": IDENTITY_POLICY,
      "mapping_inference_allowed": False,
      "model_count": len(source_models),
      "producer_repo": "alphafold3",
      "promotion_allowed": False,
      "readiness_verdict": tables.readiness_verdict,
      "residue_numbering_policy": tables.residue_numbering_policy,
      "schema_version": SCHEMA_VERSION,
      "source_models": [
          {
              "confidence_path": (
                  str(model.confidence_path) if model.confidence_path else ""
              ),
              "model_path": str(model.model_path),
              "model_rank": model.model_rank,
              "rank_index": model.rank_index,
              "ranking_score": model.ranking_score,
              "sample": model.sample,
              "seed": model.seed,
              "summary_confidence_path": (
                  str(model.summary_confidence_path)
                  if model.summary_confidence_path
                  else ""
              ),
          }
          for model in sorted(source_models, key=lambda source: source.rank_index)
      ],
      "source_output_dir": str(source_output_dir),
      "system_id": system_id,
      "table_paths": {
          "atom_identity_table": "atom_identity_table.parquet",
          "chain_identity_table": "chain_identity_table.csv",
          "confidence_geometry_correlation": "confidence_geometry_correlation.csv",
          "geometry_diversity_table": "geometry_diversity_table.csv",
          "local_class_diversity": "local_class_diversity.csv",
          "model_ranking_table": "model_ranking_table.csv",
          "pairwise_model_rmsd": "pairwise_model_rmsd.csv",
          "per_residue_variance": "per_residue_variance.csv",
          "residue_identity_table": "residue_identity_table.parquet",
      },
      "training_surface_allowed": False,
  }


def write_report(
    *,
    path: Path,
    title: str,
    manifest: Mapping[str, Any],
    metrics: Mapping[str, Any] | None = None,
    figures: Sequence[str] = (),
) -> None:
  lines = [
      f"# {title}",
      "",
      f"- System ID: `{manifest.get('system_id', '')}`",
      f"- Readiness verdict: `{manifest.get('readiness_verdict', '')}`",
      f"- Identity grade: `{manifest.get('identity_grade', '')}`",
      f"- Geometry grade: `{manifest.get('geometry_grade', '')}`",
      f"- Residue numbering policy: `{manifest.get('residue_numbering_policy', '')}`",
      f"- Atom table hash: `{manifest.get('atom_table_hash', '')}`",
      f"- Blocked reasons: `{', '.join(manifest.get('blocked_reasons', [])) or 'none'}`",
      "",
      "## Guardrails",
      "",
      "- `promotion_allowed`: false",
      "- `training_surface_allowed`: false",
      "- `contains_nmr_residuals`: false",
      "- `contains_reference_structure`: false",
      "- `mapping_inference_allowed`: false",
  ]
  if metrics:
    lines.extend(["", "## Metrics", ""])
    for key in sorted(metrics):
      lines.append(f"- `{key}`: `{metrics[key]}`")
  if figures:
    lines.extend(["", "## Figures", ""])
    for figure in figures:
      lines.append(f"![{Path(figure).stem}]({figure})")
  path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def normalize_insertion_code(value: Any) -> str:
  if value is None:
    return ""
  text = str(value)
  return "" if text in ("", ".", "?") else text


def normalize_altloc(value: Any) -> str:
  if value is None:
    return ""
  text = str(value)
  return "" if text in ("", ".", "?") else text


def canonical_atom_identity(row: Mapping[str, Any]) -> tuple[str, ...]:
  return tuple(str(row[field]) for field in CANONICAL_IDENTITY_FIELDS)


def canonical_atom_sort_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
  return (
      row["system_id"],
      row["model_rank"],
      row["chain_id"],
      _residue_number_sort_value(row["residue_number"]),
      row["insertion_code"],
      row["atom_name"],
      row["element"],
  )


def identity_without_model(row: Mapping[str, Any]) -> tuple[str, ...]:
  return (
      str(row["system_id"]),
      str(row["chain_id"]),
      str(row["residue_number"]),
      str(row["insertion_code"]),
      str(row["atom_name"]),
      str(row["element"]),
  )


def residue_identity(row: Mapping[str, Any]) -> tuple[str, ...]:
  return (
      str(row["system_id"]),
      str(row["model_rank"]),
      str(row["chain_id"]),
      str(row["residue_number"]),
      str(row["insertion_code"]),
      str(row.get("residue_name", "")),
  )


def residue_identity_without_model(row: Mapping[str, Any]) -> tuple[str, ...]:
  return (
      str(row["system_id"]),
      str(row["chain_id"]),
      str(row["residue_number"]),
      str(row["insertion_code"]),
      str(row.get("residue_name", "")),
  )


def atom_identity_hash(row: Mapping[str, Any]) -> str:
  joined = "|".join(canonical_atom_identity(row))
  return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def atom_table_hash(atom_rows: Sequence[Mapping[str, Any]]) -> str:
  payload = "\n".join(
      str(row["atom_identity_hash"])
      for row in sorted(atom_rows, key=canonical_atom_sort_key)
  )
  return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def residue_rows_from_atoms(
    atom_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  grouped: dict[tuple[str, ...], dict[str, Any]] = {}
  for row in atom_rows:
    key = residue_identity(row)
    existing = grouped.setdefault(
        key,
        {
            "system_id": row["system_id"],
            "model_rank": row["model_rank"],
            "rank_index": row["rank_index"],
            "chain_id": row["chain_id"],
            "residue_number": row["residue_number"],
            "insertion_code": row["insertion_code"],
            "residue_name": row.get("residue_name", ""),
            "label_asym_id": row.get("label_asym_id", ""),
            "label_seq_id": row.get("label_seq_id", ""),
            "atom_count": 0,
            "heavy_atom_count": 0,
        },
    )
    existing["atom_count"] += 1
    if str(row.get("element", "")).upper() != "H":
      existing["heavy_atom_count"] += 1
  return list(grouped.values())


def chain_rows_from_atoms(
    atom_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
  grouped: dict[tuple[str, str, str], dict[str, Any]] = {}
  residues: dict[tuple[str, str, str], set[tuple[str, str, str]]] = defaultdict(set)
  label_asym_ids: dict[tuple[str, str, str], set[str]] = defaultdict(set)
  for row in atom_rows:
    key = (row["system_id"], row["model_rank"], row["chain_id"])
    existing = grouped.setdefault(
        key,
        {
            "system_id": row["system_id"],
            "model_rank": row["model_rank"],
            "rank_index": row["rank_index"],
            "chain_id": row["chain_id"],
            "atom_count": 0,
            "residue_count": 0,
            "label_asym_ids": "",
        },
    )
    existing["atom_count"] += 1
    residues[key].add(
        (row["residue_number"], row["insertion_code"], row.get("residue_name", ""))
    )
    if row.get("label_asym_id"):
      label_asym_ids[key].add(str(row["label_asym_id"]))
  for key, row in grouped.items():
    row["residue_count"] = len(residues[key])
    row["label_asym_ids"] = ";".join(sorted(label_asym_ids[key]))
  return list(grouped.values())


def residue_numbering_policy(atom_rows: Sequence[Mapping[str, Any]]) -> str:
  if not atom_rows:
    return RESIDUE_NUMBERING_AUTH_PRESERVED
  per_chain: dict[tuple[str, str], list[str]] = defaultdict(list)
  seen: set[tuple[str, str, str, str]] = set()
  for row in sorted(
      atom_rows,
      key=lambda atom: (
          atom["model_rank"],
          atom["chain_id"],
          int(atom.get("atom_index", 0)),
      ),
  ):
    key = (
        row["model_rank"],
        row["chain_id"],
        row["residue_number"],
        row["insertion_code"],
    )
    if key in seen:
      continue
    seen.add(key)
    per_chain[(row["model_rank"], row["chain_id"])].append(row["residue_number"])

  found_gap = False
  for residue_numbers in per_chain.values():
    try:
      numeric = [int(value) for value in residue_numbers]
    except ValueError:
      return RESIDUE_NUMBERING_NON_MONOTONIC
    if any(curr <= prev for prev, curr in zip(numeric, numeric[1:])):
      return RESIDUE_NUMBERING_NON_MONOTONIC
    if any(curr - prev > 1 for prev, curr in zip(numeric, numeric[1:])):
      found_gap = True
  if found_gap:
    return RESIDUE_NUMBERING_CONTAINS_GAPS
  return RESIDUE_NUMBERING_AUTH_PRESERVED


def readiness_verdict(blocked_reasons: Sequence[str]) -> str:
  blocked = set(blocked_reasons)
  if not blocked:
    return READINESS_MAPPING_GRADE
  if blocked <= _REPAIRABLE_BLOCKERS:
    return READINESS_MAPPING_REPAIRABLE
  return READINESS_NOT_MAPPING_GRADE


def identity_grade(blocked_reasons: Sequence[str], policy: str) -> str:
  blocked = set(blocked_reasons)
  hard_identity_blockers = {
      BLOCKED_MISSING_ATOM_SITE_COLUMNS,
      BLOCKED_DUPLICATE_ATOM_IDENTITIES,
      BLOCKED_NONTRIVIAL_ALTLOC,
  }
  if blocked & hard_identity_blockers:
    return IDENTITY_UNSTABLE
  if policy != RESIDUE_NUMBERING_AUTH_PRESERVED:
    return IDENTITY_REPAIRABLE
  return IDENTITY_STABLE


def geometry_grade(blocked_reasons: Sequence[str]) -> str:
  blocked = set(blocked_reasons)
  if blocked & {
      BLOCKED_MISSING_ATOM_SITE_COLUMNS,
      BLOCKED_INVALID_COORDINATES,
      BLOCKED_CONFIDENCE_ATOM_COUNT_MISMATCH,
      BLOCKED_DUPLICATE_ATOM_IDENTITIES,
      BLOCKED_NONTRIVIAL_ALTLOC,
  }:
    return GEOMETRY_UNUSABLE
  if blocked:
    return GEOMETRY_PARTIALLY_COMPLETE
  return GEOMETRY_COMPLETE


def flatten_blocked_reasons(extractions: Iterable[ModelExtraction]) -> tuple[str, ...]:
  return tuple(
      sorted({reason for extraction in extractions for reason in extraction.blocked_reasons})
  )


def _atom_site_row_count(atom_site: Mapping[str, Sequence[str]]) -> int:
  lengths = {len(values) for values in atom_site.values()}
  if len(lengths) != 1:
    raise ValueError("atom_site columns have inconsistent lengths")
  return next(iter(lengths), 0)


def _value(
    table: Mapping[str, Sequence[str]], column: str, idx: int, *, default: str | None = None
) -> str:
  if column not in table:
    if default is None:
      raise KeyError(column)
    return default
  return str(table[column][idx])


def _float_or_none(value: Any) -> float | None:
  if value in ("", ".", "?", None):
    return None
  try:
    result = float(value)
  except (TypeError, ValueError):
    return None
  return result if math.isfinite(result) else None


def _load_confidence(path: Path | None) -> Mapping[str, Any]:
  if path is None or not path.exists():
    return {}
  try:
    data = read_json_maybe_zst(path)
  except (json.JSONDecodeError, OSError, zstandard.ZstdError):
    return {}
  return data if isinstance(data, Mapping) else {}


def _duplicate_identity_counts(rows: Sequence[Mapping[str, Any]]) -> tuple[int, int]:
  identities = [canonical_atom_identity(row) for row in rows]
  counts = Counter(identities)
  duplicate_count = sum(count - 1 for count in counts.values() if count > 1)
  altlocs_by_identity: dict[tuple[str, ...], set[str]] = defaultdict(set)
  for row in rows:
    altlocs_by_identity[canonical_atom_identity(row)].add(str(row.get("altloc", "")))
  duplicate_altloc_count = sum(
      1 for altlocs in altlocs_by_identity.values() if len(altlocs) > 1
  )
  return duplicate_count, duplicate_altloc_count


def _residue_number_sort_value(value: Any) -> tuple[int, Any]:
  text = str(value)
  if re.fullmatch(r"-?\d+", text):
    return (0, int(text))
  return (1, text)


def remove_path(path: Path) -> None:
  if path.is_dir():
    shutil.rmtree(path)
  elif path.exists():
    path.unlink()
