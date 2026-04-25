"""Synthetic AF3 output fixtures for NMR geometry-readiness tests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import zstandard


BACKBONE = (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))


def write_af3_job(
    root: Path,
    *,
    job_name: str = "toy",
    n_models: int = 2,
    n_residues: int = 6,
    perturb: bool = False,
    confidence_mismatch: bool = False,
    duplicate_atom: bool = False,
    altloc: str = ".",
    compressed: bool = False,
    hydrogens: bool = False,
) -> Path:
  job_dir = root / job_name
  job_dir.mkdir(parents=True, exist_ok=True)
  ranking_rows = []
  for sample in range(n_models):
    seed = 1
    score = 1.0 - sample * 0.1
    ranking_rows.append((seed, sample, score))
    sample_dir = job_dir / f"seed-{seed}_sample-{sample}"
    sample_dir.mkdir()
    prefix = f"{job_name}_seed-{seed}_sample-{sample}"
    cif_text, atom_count = make_cif(
        n_residues=n_residues,
        model_index=sample,
        perturb=perturb,
        duplicate_atom=duplicate_atom and sample == 0,
        altloc=altloc,
        hydrogens=hydrogens,
    )
    model_path = sample_dir / f"{prefix}_model.cif"
    if compressed:
      _write_zst(model_path.with_suffix(model_path.suffix + ".zst"), cif_text)
    else:
      model_path.write_text(cif_text, encoding="utf-8")
    plddt_count = atom_count - 1 if confidence_mismatch and sample == 0 else atom_count
    confidence = {
        "atom_chain_ids": ["A"] * plddt_count,
        "atom_plddts": _plddts(plddt_count, n_residues),
    }
    confidence_text = json.dumps(confidence, sort_keys=True)
    confidence_path = sample_dir / f"{prefix}_confidences.json"
    if compressed:
      _write_zst(
          confidence_path.with_suffix(confidence_path.suffix + ".zst"),
          confidence_text,
      )
    else:
      confidence_path.write_text(confidence_text, encoding="utf-8")
    (sample_dir / f"{prefix}_summary_confidences.json").write_text(
        json.dumps({"ranking_score": score}, sort_keys=True), encoding="utf-8"
    )
  with (job_dir / f"{job_name}_ranking_scores.csv").open("w", encoding="utf-8") as f:
    f.write("seed,sample,ranking_score\n")
    for seed, sample, score in ranking_rows:
      f.write(f"{seed},{sample},{score}\n")
  return job_dir


def make_cif(
    *,
    n_residues: int,
    model_index: int = 0,
    perturb: bool = False,
    duplicate_atom: bool = False,
    altloc: str = ".",
    hydrogens: bool = False,
) -> tuple[str, int]:
  rows = []
  atom_id = 1
  for res_i in range(1, n_residues + 1):
    res_name = "ALA"
    atoms: list[tuple[str, str]] = list(BACKBONE) + [("CB", "C")]
    if hydrogens:
      atoms.append(("H", "H"))
    if res_i == 3:
      res_name = "PHE"
      atoms.extend(
          [
              ("CG", "C"),
              ("CD1", "C"),
              ("CD2", "C"),
              ("CE1", "C"),
              ("CE2", "C"),
              ("CZ", "C"),
          ]
      )
    for atom_name, element in atoms:
      x = res_i * 3.0 + _atom_offset(atom_name)
      y = _atom_offset(atom_name) * 0.2
      z = 0.0
      if perturb and model_index:
        z += model_index * 0.8 if res_i > n_residues // 2 else 0.05 * model_index
        y += model_index * 0.03 * res_i
      rows.append(
          _atom_row(
              atom_id,
              element,
              atom_name,
              altloc,
              res_name,
              "A",
              res_i,
              x,
              y,
              z,
          )
      )
      atom_id += 1
  if duplicate_atom:
    rows.append(rows[0].replace("ATOM 1", f"ATOM {atom_id}", 1))
  text = "\n".join(
      [
          "data_toy",
          "loop_",
          "_atom_site.group_PDB",
          "_atom_site.id",
          "_atom_site.type_symbol",
          "_atom_site.label_atom_id",
          "_atom_site.label_alt_id",
          "_atom_site.label_comp_id",
          "_atom_site.label_asym_id",
          "_atom_site.label_entity_id",
          "_atom_site.label_seq_id",
          "_atom_site.auth_asym_id",
          "_atom_site.auth_seq_id",
          "_atom_site.pdbx_PDB_ins_code",
          "_atom_site.Cartn_x",
          "_atom_site.Cartn_y",
          "_atom_site.Cartn_z",
          "_atom_site.occupancy",
          "_atom_site.B_iso_or_equiv",
          "_atom_site.pdbx_PDB_model_num",
          *rows,
          "#",
          "",
      ]
  )
  return text, len(rows)


def make_cif_missing_required_column() -> str:
  text, _ = make_cif(n_residues=2)
  return text.replace("_atom_site.auth_seq_id\n", "")


def _atom_row(
    atom_id: int,
    element: str,
    atom_name: str,
    altloc: str,
    res_name: str,
    chain_id: str,
    res_i: int,
    x: float,
    y: float,
    z: float,
) -> str:
  return (
      f"ATOM {atom_id} {element} {atom_name} {altloc} {res_name} {chain_id} "
      f"1 {res_i} {chain_id} {res_i} ? {x:.3f} {y:.3f} {z:.3f} 1.00 80.00 1"
  )


def _atom_offset(atom_name: str) -> float:
  return sum(ord(ch) for ch in atom_name) % 7


def _plddts(count: int, n_residues: int) -> list[float]:
  values = []
  for i in range(count):
    residue = min(n_residues, i // 5 + 1)
    values.append(float(95 - residue))
  return values


def _write_zst(path: Path, text: str) -> None:
  path.write_bytes(zstandard.ZstdCompressor().compress(text.encode("utf-8")))
