#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from construct_report.pipeline import (
    build_dataset_payload,
    compute_constructs,
    write_constructs_fasta,
    write_constructs_tsv,
)
from construct_report.report import report_from_bundle

DEFAULT_PARAMS = {
    "slop": 50,
    "offset": 15,
    "minStructuredRun": 6,
    "nTerminalSnapThreshold": 15,
}

CODON_TABLE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

PDB_RESIDUE_TO_AA = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}


def clean_lines(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if line.strip()]


def parse_fasta(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                records[current_id] = "".join(parts)
            current_id = line[1:].split()[0]
            parts = []
            continue
        parts.append(line)

    if current_id:
        records[current_id] = "".join(parts)

    return records


def parse_merged_domains(text: str) -> dict[str, list[dict[str, Any]]]:
    domains: dict[str, list[dict[str, Any]]] = {}
    for line in clean_lines(text):
        protein_id, start, end, label = line.split("\t")
        domains.setdefault(protein_id, []).append(
            {
                "start": int(start),
                "end": int(end),
                "label": label,
            }
        )
    return domains


def parse_individual_domains(text: str) -> dict[str, list[dict[str, Any]]]:
    domains: dict[str, list[dict[str, Any]]] = {}
    for line in clean_lines(text):
        fields = line.split("\t")
        if len(fields) == 4:
            protein_id, start, end, label = fields
            score = None
        elif len(fields) >= 6:
            protein_id, start, end, label, _pfam_id, score = fields[:6]
        else:
            raise ValueError("Individual domain file must have either 4 or at least 6 tab-separated columns.")
        domains.setdefault(protein_id, []).append(
            {
                "start": int(start),
                "end": int(end),
                "label": label,
                "source": "domains.individual",
                "score": float(score) if score not in (None, "") else None,
            }
        )
    return domains


def parse_custom_ranges(text: str) -> dict[str, list[dict[str, Any]]]:
    lines = clean_lines(text)
    if not lines:
        return {}

    id_fields = {"gene", "id", "protein_id", "protein", "sequence_id", "seqid"}
    range_fields = {"range", "picked_range", "construct_range", "aa_range", "custom_range"}
    start_fields = {"start", "aa_start"}
    end_fields = {"end", "aa_end"}
    label_fields = {"label", "name", "source", "recommended", "range_name", "picked_range_name"}
    missing_values = {"", "na", "n/a", "none", "null", "."}

    def is_missing(value: str) -> bool:
        return value.strip().lower() in missing_values

    def parse_range_text(value: str, line_number: int) -> tuple[int, int]:
        match = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", value)
        if not match:
            raise ValueError(f"Invalid range {value!r} on line {line_number}. Expected START-END.")
        start, end = int(match.group(1)), int(match.group(2))
        if end < start:
            raise ValueError(f"Invalid range {value!r} on line {line_number}. End precedes start.")
        return start, end

    def add_entry(
        rows: dict[str, list[dict[str, Any]]],
        protein_id: str,
        start: int,
        end: int,
        label: str | None,
    ) -> None:
        rows.setdefault(protein_id, []).append(
            {
                "start": start,
                "end": end,
                "label": label or "custom range",
            }
        )

    header_fields = [field.strip().lower() for field in lines[0].split("\t")]
    has_header = any(
        field in id_fields | range_fields | start_fields | end_fields | label_fields
        for field in header_fields
    )
    rows: dict[str, list[dict[str, Any]]] = {}

    if has_header:
        id_index = next((i for i, field in enumerate(header_fields) if field in id_fields), None)
        range_index = next((i for i, field in enumerate(header_fields) if field in range_fields), None)
        start_index = next((i for i, field in enumerate(header_fields) if field in start_fields), None)
        end_index = next((i for i, field in enumerate(header_fields) if field in end_fields), None)
        label_index = next((i for i, field in enumerate(header_fields) if field in label_fields), None)

        if id_index is None:
            raise ValueError("Custom ranges header is missing an ID column.")
        if range_index is None and (start_index is None or end_index is None):
            raise ValueError(
                "Custom ranges header must contain either a range column or both start/end columns."
            )

        for line_number, line in enumerate(lines[1:], start=2):
            fields = line.split("\t")
            protein_id = fields[id_index].strip() if id_index < len(fields) else ""
            if not protein_id:
                continue

            label = (
                fields[label_index].strip()
                if label_index is not None and label_index < len(fields)
                else None
            )
            range_value = fields[range_index].strip() if range_index is not None and range_index < len(fields) else ""
            start_value = fields[start_index].strip() if start_index is not None and start_index < len(fields) else ""
            end_value = fields[end_index].strip() if end_index is not None and end_index < len(fields) else ""

            if not is_missing(range_value):
                start, end = parse_range_text(range_value, line_number)
            else:
                if is_missing(start_value) and is_missing(end_value):
                    continue
                if is_missing(start_value) or is_missing(end_value):
                    raise ValueError(f"Custom ranges line {line_number} is missing start/end values.")
                start = int(start_value)
                end = int(end_value)
                if end < start:
                    raise ValueError(f"Custom ranges line {line_number} has end before start.")
            add_entry(rows, protein_id, start, end, label)

        return rows

    for line_number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(
                f"Custom ranges line {line_number} must have at least 2 tab-separated fields."
            )

        protein_id = fields[0].strip()
        if not protein_id:
            continue

        if len(fields) >= 3 and fields[1].strip().isdigit() and fields[2].strip().isdigit():
            start = int(fields[1].strip())
            end = int(fields[2].strip())
            if end < start:
                raise ValueError(f"Custom ranges line {line_number} has end before start.")
            label = fields[3].strip() if len(fields) >= 4 and fields[3].strip() else None
        else:
            start, end = parse_range_text(fields[1].strip(), line_number)
            label = fields[2].strip() if len(fields) >= 3 and fields[2].strip() else None

        add_entry(rows, protein_id, start, end, label)

    return rows


def parse_metadata_table(
    text: str,
    allowed_ids: set[str] | None = None,
) -> dict[str, Any]:
    if not text.strip():
        return {"columns": [], "rows": [], "idKey": None}

    sample = "\n".join(text.splitlines()[:5])
    tab_count = sample.count("\t")
    comma_count = sample.count(",")
    delimiter = "\t" if tab_count >= comma_count else ","
    reader = csv.reader(io.StringIO(text), delimiter=delimiter)
    rows = [[cell.strip() for cell in row] for row in reader if any(cell.strip() for cell in row)]
    if not rows:
        return {"columns": [], "rows": [], "idKey": None}

    raw_header = rows[0]
    width = max(len(raw_header), *(len(row) for row in rows[1:]))
    padded_header = raw_header + [""] * (width - len(raw_header))

    header_names: list[str] = []
    seen_names: dict[str, int] = {}
    for index, raw_name in enumerate(padded_header):
        base_name = raw_name.strip().lstrip("\ufeff") or f"column_{index + 1}"
        occurrence = seen_names.get(base_name, 0)
        seen_names[base_name] = occurrence + 1
        header_names.append(base_name if occurrence == 0 else f"{base_name}_{occurrence + 1}")

    lowered = [name.lower() for name in header_names]
    id_index = lowered.index("gene") if "gene" in lowered else (2 if width >= 3 else 0)
    ordered_indices = [id_index] + [index for index in range(width) if index != id_index]
    ordered_columns = [
        {
            "key": header_names[index],
            "label": header_names[index],
        }
        for index in ordered_indices
    ]
    id_key = header_names[id_index]

    table_rows: list[dict[str, str]] = []
    for raw_row in rows[1:]:
        padded_row = raw_row + [""] * (width - len(raw_row))
        row_map = {header_names[index]: padded_row[index] for index in range(width)}
        protein_id = row_map.get(id_key, "").strip()
        if not protein_id:
            continue
        if allowed_ids is not None and protein_id not in allowed_ids:
            continue
        table_rows.append({column["key"]: row_map.get(column["key"], "") for column in ordered_columns})

    return {
        "columns": ordered_columns,
        "rows": table_rows,
        "idKey": id_key,
    }


def locate_sequence(target: str, query: str) -> int | None:
    if target == query:
        return 1
    match_index = target.find(query)
    return match_index + 1 if match_index >= 0 else None


def parse_conservation_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    if not lines:
        raise ValueError("Conservation file is empty.")
    if len(lines) < 2:
        raise ValueError(f"Conservation file has only 1 line (need sequence + scores): {lines[0][:60]!r}")

    sequence = lines[0]
    if len(lines) == 2:
        try:
            values = [float(value) for value in lines[1].split(",")]
        except ValueError as exc:
            raise ValueError(f"Conservation scores line is not comma-separated floats: {exc}") from exc
        return {"sequence": sequence, "values": values}

    alignment = list(lines[1])
    raw_scores = [float(value) for value in lines[2].split(",")]
    values = [score for score, residue in zip(raw_scores, alignment) if residue != "-"]
    return {"sequence": sequence, "values": values}


def parse_iupred_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    if not lines:
        raise ValueError("IUPred file is empty.")

    fields = lines[0].split("\t")
    if len(fields) < 2:
        raise ValueError("IUPred file must contain sequence and comma-separated scores.")

    sequence = fields[0].strip()
    if not sequence:
        raise ValueError("IUPred file is missing the sequence column.")

    try:
        values = [float(value) for value in fields[1].split(",") if value.strip()]
    except ValueError as exc:
        raise ValueError(f"IUPred scores are not comma-separated floats: {exc}") from exc

    if len(values) != len(sequence):
        raise ValueError(
            f"IUPred score count ({len(values)}) does not match sequence length ({len(sequence)})."
        )

    return {
        "sequence": sequence,
        "values": values,
        "maxValue": 1.0,
    }


def parse_structure_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    sequence_index = next((i for i, line in enumerate(lines) if line.endswith("@sequence")), -1)
    structure_index = next((i for i, line in enumerate(lines) if line.endswith("@secondary")), -1)

    if sequence_index < 0 or structure_index < 0:
        raise ValueError("Structure file is missing sequence or secondary blocks.")

    sequence = lines[sequence_index + 1] if sequence_index + 1 < len(lines) else ""
    secondary = lines[structure_index + 1] if structure_index + 1 < len(lines) else ""
    return {"sequence": sequence, "values": list(secondary)}


def parse_pdb_model(text: str) -> dict[str, Any]:
    residue_numbers: list[int] = []
    sequence: list[str] = []
    seen: set[tuple[str, str, str]] = set()

    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue

        chain_id = line[21]
        residue_number_text = line[22:26].strip()
        insertion_code = line[26]
        residue_key = (chain_id, residue_number_text, insertion_code)
        if residue_key in seen:
            continue

        try:
            residue_number = int(residue_number_text)
        except ValueError:
            continue

        seen.add(residue_key)
        residue_numbers.append(residue_number)
        sequence.append(PDB_RESIDUE_TO_AA.get(line[17:20].strip().upper(), "X"))

    return {
        "residueNumbers": residue_numbers,
        "sequence": "".join(sequence),
    }


def parse_plddt_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    if not lines:
        raise ValueError("pLDDT file is empty.")
    header = lines[0].split("\t")
    if "Positions" not in header or "rank_0" not in header:
        raise ValueError("pLDDT TSV must have 'Positions' and 'rank_0' columns.")
    pos_col = header.index("Positions")
    val_col = header.index("rank_0")
    positions: list[int] = []
    values: list[float | None] = []
    for line in lines[1:]:
        fields = line.split("\t")
        try:
            positions.append(int(fields[pos_col]))
        except (ValueError, IndexError):
            continue
        try:
            raw = fields[val_col].strip()
            values.append(float(raw) if raw else None)
        except (ValueError, IndexError):
            values.append(None)
    return {"positions": positions, "values": values}


def map_plddt_track(protein_sequence: str, protein_start: int, raw: dict[str, Any]) -> dict[str, Any]:
    full_values: list[float | None] = [None] * len(protein_sequence)
    for pos, val in zip(raw["positions"], raw["values"]):
        idx = protein_start - 1 + (pos - 1)
        if 0 <= idx < len(full_values):
            full_values[idx] = val
    non_null = sum(1 for v in full_values if v is not None)
    return {
        "type": "numeric",
        "label": "pLDDT",
        "offset": protein_start,
        "values": full_values,
        "coverage": non_null / len(protein_sequence) if protein_sequence else 0,
        "compatible": non_null > 0,
        "maxValue": 100.0,
    }


def map_numeric_track(protein_sequence: str, label: str, raw_track: dict[str, Any]) -> dict[str, Any]:
    offset = locate_sequence(protein_sequence, raw_track["sequence"])
    if not offset:
        return {
            "type": "numeric",
            "label": label,
            "sequence": raw_track["sequence"],
            "offset": None,
            "values": [None for _ in protein_sequence],
            "coverage": 0,
            "compatible": False,
            "reason": "No exact in-sequence placement found.",
        }

    values: list[float | None] = [None for _ in protein_sequence]
    for index, value in enumerate(raw_track["values"]):
        target_index = offset - 1 + index
        if 0 <= target_index < len(values):
            values[target_index] = float(value)

    mapped_track = {
        "type": "numeric",
        "label": label,
        "sequence": raw_track["sequence"],
        "offset": offset,
        "values": values,
        "coverage": len(raw_track["values"]) / len(protein_sequence),
        "compatible": True,
    }
    if raw_track.get("maxValue") is not None:
        mapped_track["maxValue"] = raw_track["maxValue"]
    return mapped_track


def map_structure_track(protein_sequence: str, label: str, raw_track: dict[str, Any]) -> dict[str, Any]:
    offset = locate_sequence(protein_sequence, raw_track["sequence"])
    if not offset:
        return {
            "type": "structure",
            "label": label,
            "sequence": raw_track["sequence"],
            "offset": None,
            "values": [None for _ in protein_sequence],
            "coverage": 0,
            "compatible": False,
            "reason": "No exact in-sequence placement found.",
        }

    values: list[str | None] = [None for _ in protein_sequence]
    for index, value in enumerate(raw_track["values"]):
        target_index = offset - 1 + index
        if 0 <= target_index < len(values):
            values[target_index] = str(value)

    return {
        "type": "structure",
        "label": label,
        "sequence": raw_track["sequence"],
        "offset": offset,
        "values": values,
        "coverage": len(raw_track["values"]) / len(protein_sequence),
        "compatible": True,
    }


def build_structure_model_mapping(
    protein_sequence: str,
    model_text: str,
    local_track: dict[str, Any] | None,
    fallback_track: dict[str, Any] | None,
) -> dict[str, Any]:
    parsed_model = parse_pdb_model(model_text)
    residue_numbers = parsed_model["residueNumbers"]
    model_sequence = parsed_model["sequence"]
    residue_count = len(residue_numbers)

    protein_start = None
    mapping_source = None
    model_length = residue_count

    # Prefer the actual PDB sequence placement when it can be matched.
    # DSSP/UniProt structure tracks often cover only the structured subset of a
    # full-length model, and using their offset to anchor the whole PDB shifts
    # the displayed model range incorrectly.
    if model_sequence:
        sequence_offset = locate_sequence(protein_sequence, model_sequence)
        if sequence_offset:
            protein_start = sequence_offset
            mapping_source = "PDB sequence"
            model_length = min(residue_count, len(model_sequence))
    if protein_start is None:
        for candidate in (local_track, fallback_track):
            if candidate and candidate.get("compatible") and candidate.get("offset"):
                protein_start = int(candidate["offset"])
                mapping_source = candidate["label"]
                model_length = min(
                    residue_count,
                    len(candidate.get("sequence", "")) or residue_count,
                )
                break
    if protein_start is None and residue_count == len(protein_sequence):
        protein_start = 1
        mapping_source = "full-length fallback"
        model_length = residue_count

    protein_end = protein_start + model_length - 1 if protein_start and model_length else None

    return {
        "mappingSource": mapping_source,
        "proteinStart": protein_start,
        "proteinEnd": protein_end,
        "modelLength": model_length,
        "residueCount": residue_count,
        "residueStart": residue_numbers[0] if residue_numbers else None,
        "residueEnd": residue_numbers[model_length - 1] if residue_numbers and model_length else None,
        "residueNumbers": residue_numbers[:model_length] if model_length else [],
    }


_CRYST1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"


def compute_dssp_from_pdb(pdb_path: Path, pdb_text: str) -> dict[str, Any] | None:
    """Run DSSP on an AlphaFold PDB and return a raw_track dict for map_structure_track."""
    try:
        import tempfile
        import os
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP

        fixed = _CRYST1 + pdb_text
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
            tmp.write(fixed)
            tmp_path = tmp.name
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("model", tmp_path)
            model = structure[0]
            dssp = DSSP(model, tmp_path)
        finally:
            os.unlink(tmp_path)

        seq, ss = [], []
        for key in dssp.keys():
            d = dssp[key]
            aa = d[1] if d[1] != "X" else "?"
            s = d[2] if d[2] != " " else "-"
            seq.append(aa)
            ss.append(s)

        if not seq:
            return None
        return {"sequence": "".join(seq), "values": ss}
    except Exception as exc:
        print(f"WARNING: DSSP computation failed for {pdb_path.name}: {exc}", file=sys.stderr)
        return None


def dssp_matches_pdb_model(raw_track: dict[str, Any] | None, pdb_text: str) -> bool:
    if not raw_track:
        return False

    dssp_sequence = str(raw_track.get("sequence") or "")
    if not dssp_sequence:
        return False

    parsed_model = parse_pdb_model(pdb_text)
    model_sequence = parsed_model.get("sequence") or ""
    if not model_sequence:
        return False

    return dssp_sequence == model_sequence


def build_evidence_for_protein(
    protein_id: str,
    protein_sequence: str,
    evidence_files: list[Path],
    structure_files: list[Path],
    plddt_files: list[Path] | None = None,
) -> dict[str, Any]:
    bundle: dict[str, Any] = {}
    raw_structure_dssp: dict[str, Any] | None = None

    # Evidence support is intentionally explicit for now: only these known
    # subdirectory names are loaded into report tracks.
    for path in evidence_files:
        if not path.name.startswith(protein_id):
            continue

        content = path.read_text(encoding="utf-8")
        kind = path.parent.name

        try:
            if kind == "conservation":
                bundle["conservation"] = map_numeric_track(
                    protein_sequence,
                    "Domain conservation",
                    parse_conservation_file(content),
                )
            elif kind == "conservation_full":
                bundle["conservationFull"] = map_numeric_track(
                    protein_sequence,
                    "Full-length conservation",
                    parse_conservation_file(content),
                )
            elif kind == "iupred":
                bundle["iupred"] = map_numeric_track(
                    protein_sequence,
                    "IUPred disorder",
                    parse_iupred_file(content),
                )
            elif kind == "structure_uniprot":
                bundle["structureUniprot"] = map_structure_track(
                    protein_sequence,
                    "UniProt structure",
                    parse_structure_file(content),
                )
            elif kind == "structure_dssp":
                raw_structure_dssp = parse_structure_file(content)
                bundle["structureDssp"] = map_structure_track(
                    protein_sequence,
                    "Local DSSP structure",
                    raw_structure_dssp,
                )
        except (ValueError, IndexError) as exc:
            print(f"WARNING: skipping {path}: {exc}", file=sys.stderr)

    matching_structures = sorted(
        (path for path in structure_files if path.name.startswith(protein_id)),
        key=lambda path: (path.stem != protein_id, path.name),
    )
    if matching_structures:
        structure_path = matching_structures[0]
        structure_text = structure_path.read_text(encoding="utf-8")

        should_recompute_dssp = False
        if "structureDssp" not in bundle:
            should_recompute_dssp = True
        elif not dssp_matches_pdb_model(raw_structure_dssp, structure_text):
            print(
                (
                    f"WARNING: precomputed DSSP for {protein_id} does not match "
                    f"{structure_path.name}; recomputing from PDB"
                ),
                file=sys.stderr,
            )
            should_recompute_dssp = True

        if should_recompute_dssp:
            raw_dssp = compute_dssp_from_pdb(structure_path, structure_text)
            if raw_dssp:
                bundle["structureDssp"] = map_structure_track(
                    protein_sequence, "Local DSSP structure", raw_dssp
                )

        bundle["structureModel"] = {
            "format": structure_path.suffix.lstrip(".") or "pdb",
            "source": structure_path.name,
            "text": structure_text,
            "mapping": build_structure_model_mapping(
                protein_sequence,
                structure_text,
                bundle.get("structureDssp"),
                bundle.get("structureUniprot"),
            ),
        }

    if plddt_files:
        matching_plddt = sorted(
            (
                path
                for path in plddt_files
                if path.name.startswith(protein_id) and path.parent.name == "structures"
            ),
            key=lambda path: (path.stem != protein_id, path.name),
        )
        if matching_plddt:
            try:
                raw_plddt = parse_plddt_file(matching_plddt[0].read_text(encoding="utf-8"))
                protein_start = (
                    bundle.get("structureModel", {}).get("mapping", {}).get("proteinStart") or 1
                )
                bundle["plddt"] = map_plddt_track(protein_sequence, protein_start, raw_plddt)
            except (ValueError, KeyError):
                pass

    return bundle


def load_dataset_bundle(
    pep_path: Path,
    cds_path: Path,
    domains_path: Path | None,
    domains_individual_path: Path | None = None,
    evidence_root: Path | None = None,
) -> dict[str, Any]:
    proteins = parse_fasta(pep_path.read_text(encoding="utf-8"))
    cds = parse_fasta(cds_path.read_text(encoding="utf-8"))
    merged_domains = (
        parse_merged_domains(domains_path.read_text(encoding="utf-8"))
        if domains_path and domains_path.exists()
        else {}
    )
    individual_domains = (
        parse_individual_domains(domains_individual_path.read_text(encoding="utf-8"))
        if domains_individual_path and domains_individual_path.exists()
        else {}
    )
    evidence_files = (
        sorted(evidence_root.rglob("*.out"))
        if evidence_root and evidence_root.exists()
        else []
    )
    structure_files = (
        sorted(evidence_root.rglob("*.pdb"))
        if evidence_root and evidence_root.exists()
        else []
    )
    plddt_files = (
        sorted(evidence_root.rglob("*.tsv"))
        if evidence_root and evidence_root.exists()
        else []
    )
    protein_ids = set(proteins)
    cds_ids = set(cds)
    shared_ids = protein_ids & cds_ids
    allowed_protein_ids = set(individual_domains) if individual_domains else None
    evidence_track_ids = {
        protein_id
        for protein_id in protein_ids
        if any(path.name.startswith(protein_id) for path in evidence_files)
    }
    structure_model_ids = {
        protein_id
        for protein_id in protein_ids
        if any(path.name.startswith(protein_id) for path in structure_files)
    }

    entries: list[dict[str, Any]] = []
    for protein_id in sorted(proteins):
        if allowed_protein_ids is not None and protein_id not in allowed_protein_ids:
            continue

        cds_sequence = cds.get(protein_id)
        if not cds_sequence:
            continue

        protein_sequence = proteins[protein_id]
        entries.append(
            {
                "id": protein_id,
                "proteinSequence": protein_sequence,
                "cdsSequence": cds_sequence,
                "mergedDomains": merged_domains.get(protein_id, []),
                "individualDomains": individual_domains.get(protein_id, []),
                "evidence": build_evidence_for_protein(
                    protein_id,
                    protein_sequence,
                    evidence_files,
                    structure_files,
                    plddt_files,
                ),
            }
        )

    summary = {
        "pepRecords": len(proteins),
        "cdsRecords": len(cds),
        "sharedPepCds": len(shared_ids),
        "pepWithoutCds": len(protein_ids - cds_ids),
        "cdsWithoutPep": len(cds_ids - protein_ids),
        "mergedDomainProteins": len(merged_domains),
        "individualDomainProteins": len(individual_domains),
        "filterRequiresIndividualDomains": allowed_protein_ids is not None,
        "excludedByIndividualDomains": len(shared_ids - allowed_protein_ids) if allowed_protein_ids is not None else 0,
        "evidenceTrackFiles": len(evidence_files),
        "structureModelFiles": len(structure_files),
        "evidenceTrackProteins": len(evidence_track_ids),
        "structureModelProteins": len(structure_model_ids),
        "evidenceProteinsAny": len(evidence_track_ids | structure_model_ids),
        "keptProteins": len(entries),
        "keptWithAnyEvidence": sum(1 for entry in entries if entry["evidence"]),
        "keptWithTrackEvidence": sum(
            1
            for entry in entries
            if any(key != "structureModel" for key in entry["evidence"])
        ),
        "keptWithStructureModels": sum(
            1 for entry in entries if entry["evidence"].get("structureModel")
        ),
    }

    return {"entries": entries, "summary": summary}


def load_dataset(
    pep_path: Path,
    cds_path: Path,
    domains_path: Path | None,
    domains_individual_path: Path | None = None,
    evidence_root: Path | None = None,
) -> list[dict[str, Any]]:
    return load_dataset_bundle(
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        domains_individual_path=domains_individual_path,
        evidence_root=evidence_root,
    )["entries"]


def resolve_path(base_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def format_path_for_display(path: Path, base_dir: Path) -> str:
    try:
        return str(path.relative_to(base_dir))
    except ValueError:
        return str(path)


def summarize_dataset_for_display(summary: dict[str, Any]) -> str:
    parts = [
        f"pep {summary['pepRecords']}",
        f"cds {summary['cdsRecords']}",
        f"shared {summary['sharedPepCds']}",
        f"kept {summary['keptProteins']}",
    ]
    if summary["evidenceProteinsAny"] or summary["evidenceTrackFiles"] or summary["structureModelFiles"]:
        parts.append(f"evidence {summary['evidenceProteinsAny']}")
    if summary["keptWithStructureModels"]:
        parts.append(f"structures {summary['keptWithStructureModels']}")
    return " | ".join(parts)


def format_run_summary_lines(summary: dict[str, Any]) -> list[str]:
    lines = [
        "Dataset summary:",
        f"  pep records: {summary['pepRecords']}",
        f"  cds records: {summary['cdsRecords']}",
        f"  pep+cds overlap: {summary['sharedPepCds']}",
        f"  pep without cds: {summary['pepWithoutCds']}",
        f"  cds without pep: {summary['cdsWithoutPep']}",
        f"  merged-domain proteins: {summary['mergedDomainProteins']}",
        f"  individual-domain proteins: {summary['individualDomainProteins']}",
        (
            "  individual-domain filter: on"
            if summary["filterRequiresIndividualDomains"]
            else "  individual-domain filter: off"
        ),
        f"  excluded by individual domains: {summary['excludedByIndividualDomains']}",
        f"  evidence track files: {summary['evidenceTrackFiles']}",
        f"  structure model files: {summary['structureModelFiles']}",
        f"  proteins referenced by evidence: {summary['evidenceProteinsAny']}",
        f"  proteins kept: {summary['keptProteins']}",
        f"  kept with any evidence: {summary['keptWithAnyEvidence']}",
        f"  kept with track evidence: {summary['keptWithTrackEvidence']}",
        f"  kept with structure models: {summary['keptWithStructureModels']}",
    ]
    return lines


def build_params(
    slop: int,
    offset: int,
    min_structured_run: int,
    n_terminal_snap_threshold: int,
) -> dict[str, int]:
    return {
        "slop": max(1, slop),
        "offset": max(1, offset),
        "minStructuredRun": max(1, min_structured_run),
        "nTerminalSnapThreshold": max(1, n_terminal_snap_threshold),
    }


def resolve_cli_inputs(
    *,
    cwd: Path,
    project_root: Path,
    input_dir_arg: str | None,
    pep_arg: str | None,
    cds_arg: str | None,
    domains_arg: str | None,
    domains_individual_arg: str | None,
    custom_ranges_arg: str | None = None,
    metadata_arg: str | None = None,
    evidence_dir_arg: str | None,
) -> dict[str, Any]:
    input_dir = (
        resolve_path(cwd, input_dir_arg)
        if input_dir_arg
        else (project_root / "examples").resolve()
    )
    primary_explicit = any(value is not None for value in (pep_arg, cds_arg, domains_arg))
    pep_path = resolve_path(cwd, pep_arg) or (input_dir / "proteins.fasta")
    cds_path = resolve_path(cwd, cds_arg) or (input_dir / "cds.fasta")

    if domains_arg:
        domains_path = resolve_path(cwd, domains_arg)
    elif primary_explicit:
        domains_path = None
    else:
        default_domains = input_dir / "domains.bed"
        domains_path = default_domains if default_domains.exists() else None

    domains_individual_path = (
        resolve_path(cwd, domains_individual_arg)
        if domains_individual_arg
        else (
            None
            if primary_explicit
            else (
                (input_dir / "domains.individual.bed")
                if (input_dir / "domains.individual.bed").exists()
                else (input_dir / "domains.individual.tab")
            )
        )
    )
    custom_ranges_path = (
        resolve_path(cwd, custom_ranges_arg)
        if custom_ranges_arg
        else (
            None
            if primary_explicit
            else (
                (input_dir / "custom_ranges.tsv")
                if (input_dir / "custom_ranges.tsv").exists()
                else ((input_dir / "custom_ranges.tab") if (input_dir / "custom_ranges.tab").exists() else None)
            )
        )
    )
    metadata_path = (
        resolve_path(cwd, metadata_arg)
        if metadata_arg
        else (
            None
            if primary_explicit
            else (
                (input_dir / "metadata.tsv")
                if (input_dir / "metadata.tsv").exists()
                else (
                    (input_dir / "metadata.csv")
                    if (input_dir / "metadata.csv").exists()
                    else (
                        (input_dir / "TF_list.csv")
                        if (input_dir / "TF_list.csv").exists()
                        else (
                            (project_root / "metadata.tsv")
                            if (project_root / "metadata.tsv").exists()
                            else (
                                (project_root / "metadata.csv")
                                if (project_root / "metadata.csv").exists()
                                else ((project_root / "TF_list.csv") if (project_root / "TF_list.csv").exists() else None)
                            )
                        )
                    )
                )
            )
        )
    )
    evidence_root = (
        resolve_path(cwd, evidence_dir_arg)
        if evidence_dir_arg
        else (None if primary_explicit else input_dir / "evidence")
    )
    explicit_sources = any(
        value is not None
        for value in (
            pep_arg,
            cds_arg,
            domains_arg,
            domains_individual_arg,
            metadata_arg,
            evidence_dir_arg,
        )
    )

    return {
        "input_dir": input_dir,
        "pep_path": pep_path,
        "cds_path": cds_path,
        "domains_path": domains_path,
        "domains_individual_path": domains_individual_path,
        "custom_ranges_path": custom_ranges_path,
        "metadata_path": metadata_path,
        "evidence_root": evidence_root,
        "explicit_sources": explicit_sources,
    }


def build_input_summary(
    *,
    cwd: Path,
    input_dir: Path,
    pep_path: Path,
    cds_path: Path,
    domains_path: Path | None,
    explicit_sources: bool,
) -> str:
    if not explicit_sources:
        return format_path_for_display(input_dir, cwd)

    return " | ".join(
        [
            f"pep={format_path_for_display(pep_path, cwd)}",
            f"cds={format_path_for_display(cds_path, cwd)}",
            f"domains={format_path_for_display(domains_path, cwd) if domains_path else 'none'}",
        ]
    )


def _pluralize(count: int, noun: str) -> str:
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def format_evidence_discovery_lines(
    *,
    evidence_root: Path | None,
) -> list[str]:
    if not evidence_root or not evidence_root.exists():
        return []

    def count_matching(subdir_name: str, pattern: str) -> int:
        subdir = evidence_root / subdir_name
        return len(list(subdir.glob(pattern))) if subdir.exists() else 0

    lines = ["  evidence found:"]
    out_kinds = [
        ("conservation", "conservation"),
        ("conservation_full", "conservation_full"),
        ("iupred", "iupred"),
        ("structure_dssp", "structure_dssp"),
        ("structure_uniprot", "structure_uniprot"),
    ]
    for subdir_name, label in out_kinds:
        count = count_matching(subdir_name, "*.out")
        detail = f"{_pluralize(count, '.out file')}" if count else "absent"
        lines.append(f"    {label}: {detail}")

    structure_dir = evidence_root / "structures"
    if structure_dir.exists():
        pdb_count = len(list(structure_dir.glob("*.pdb")))
        plddt_count = len(list(structure_dir.glob("*.tsv")))
        if pdb_count or plddt_count:
            parts = []
            if pdb_count:
                parts.append(_pluralize(pdb_count, "PDB"))
            if plddt_count:
                parts.append(_pluralize(plddt_count, "pLDDT TSV"))
            lines.append(f"    structures: {', '.join(parts)}")
        else:
            lines.append("    structures: present, but no .pdb or .tsv files found")
    else:
        lines.append("    structures: absent")

    return lines


def format_resolved_input_lines(
    *,
    cwd: Path,
    paths: dict[str, Any],
) -> list[str]:
    heading = (
        "Resolved inputs:"
        if paths["explicit_sources"]
        else "Picked up from input directory:"
    )
    lines = [
        heading,
        f"  input dir: {format_path_for_display(paths['input_dir'], cwd)}",
        f"  proteins: {format_path_for_display(paths['pep_path'], cwd)}",
        f"  cds: {format_path_for_display(paths['cds_path'], cwd)}",
    ]

    optional_items = [
        ("domains", paths["domains_path"]),
        ("individual domains", paths["domains_individual_path"]),
        ("custom ranges", paths["custom_ranges_path"]),
        ("metadata", paths["metadata_path"]),
        ("evidence dir", paths["evidence_root"]),
    ]
    for label, path in optional_items:
        if path and path.exists():
            lines.append(f"  {label}: {format_path_for_display(path, cwd)}")
            if label == "evidence dir":
                lines.extend(
                    format_evidence_discovery_lines(
                        evidence_root=path,
                    )
                )

    return lines


def append_custom_ranges_to_payload(
    payload: dict[str, Any],
    custom_ranges_path: Path | None,
    cwd: Path,
) -> dict[str, Any]:
    if not custom_ranges_path:
        return payload

    payload["customRanges"] = parse_custom_ranges(
        custom_ranges_path.read_text(encoding="utf-8")
    )
    custom_display = format_path_for_display(custom_ranges_path, cwd)
    current_summary = str(payload.get("inputSummary", "") or "")
    payload["inputSummary"] = (
        f"{current_summary} | custom={custom_display}"
        if current_summary
        else f"custom={custom_display}"
    )
    return payload


def append_metadata_to_payload(
    payload: dict[str, Any],
    metadata_path: Path | None,
    cwd: Path,
) -> dict[str, Any]:
    if not metadata_path:
        return payload

    dataset_ids = {
        str(entry.get("id"))
        for entry in payload.get("dataset", [])
        if entry.get("id")
    }
    metadata_table = parse_metadata_table(
        metadata_path.read_text(encoding="utf-8"),
        allowed_ids=dataset_ids or None,
    )
    if metadata_table.get("columns") and metadata_table.get("rows"):
        metadata_table["source"] = format_path_for_display(metadata_path, cwd)
        payload["metadataTable"] = metadata_table
        current_summary = str(payload.get("inputSummary", "") or "")
        metadata_display = format_path_for_display(metadata_path, cwd)
        payload["inputSummary"] = (
            f"{current_summary} | metadata={metadata_display}"
            if current_summary
            else f"metadata={metadata_display}"
        )
    return payload


def render_html(payload_json: str) -> str:
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Construct Design HTML Report</title>
  <style>
    :root {
      --bg: #e4e4e4;
      --panel: #ffffff;
      --ink: #1a1a1a;
      --muted: #666666;
      --border: #aaaaaa;
      --border-light: #cccccc;
      --border-inner: #d8d8d8;
      --accent: #2c3e50;
      --good: #27ae60;
      --bad: #c0392b;
      --amber: #d4943a;
      --info: #2868a0;
      --shadow: 1px 1px 2px rgba(0, 0, 0, 0.1);
      --font-ui: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      --font-mono: "SFMono-Regular", "Menlo", "Consolas", "Courier New", monospace;
    }

    * {
      box-sizing: border-box;
    }

    body {
      margin: 0;
      min-width: 320px;
      background: var(--bg);
      color: var(--ink);
      font-family: var(--font-ui);
      font-size: 15px;
      line-height: 1.5;
    }

    a {
      color: inherit;
      text-decoration: none;
    }

    button,
    input,
    select {
      font: inherit;
    }

    .app-shell {
      padding: 12px;
      min-height: 100vh;
    }

    .app-toolbar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      padding: 6px 10px;
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 3px;
      box-shadow: var(--shadow);
      margin-bottom: 8px;
    }

    .app-toolbar-title {
      display: flex;
      flex-direction: column;
      gap: 1px;
    }

    .app-toolbar-title strong {
      font-size: 16px;
    }

    .app-toolbar-title span {
      font-size: 12px;
      color: var(--muted);
    }

    .app-toolbar-actions {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
    }

    .app-tabs {
      display: flex;
      gap: 6px;
      margin-bottom: 8px;
    }

    .app-tab {
      appearance: none;
      border: 1px solid var(--border);
      background: #f2f2f2;
      color: #444;
      font: inherit;
      font-size: 11px;
      line-height: 1;
      padding: 7px 10px;
      border-radius: 999px;
      cursor: pointer;
    }

    .app-tab:hover {
      background: #eceff2;
    }

    .app-tab-active {
      background: #dde7ef;
      border-color: var(--accent);
      color: #18384e;
      font-weight: 700;
    }

    .app-page-hidden {
      display: none !important;
    }

    .action-button {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-height: 26px;
      padding: 0 10px;
      border: 1px solid var(--border);
      background: #f5f5f5;
      border-radius: 2px;
      cursor: pointer;
    }

    .action-button:hover {
      background: #ebf1f6;
    }

    .action-button-secondary {
      background: #eef4f8;
    }

    .action-button-muted {
      color: var(--muted);
      background: #f0f0f0;
      cursor: default;
    }

    .workspace {
      display: grid;
      grid-template-columns: minmax(240px, 290px) minmax(0, 1fr);
      gap: 8px;
    }

    .panel {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 3px;
      box-shadow: var(--shadow);
    }

    .batch-panel {
      padding: 0;
      display: flex;
      flex-direction: column;
      min-height: calc(100vh - 80px);
      overflow: hidden;
    }

    .batch-panel-header {
      padding: 8px 10px 7px;
      border-bottom: 1px solid var(--border-light);
      background: #f4f4f4;
    }

    .batch-panel-header h2 {
      margin: 0;
      font-size: 12px;
    }

    .batch-panel-header p {
      margin: 3px 0 0;
      font-size: 10px;
      color: var(--muted);
    }

    .batch-toolbar {
      display: flex;
      gap: 4px;
      padding: 5px 6px;
      border-bottom: 1px solid var(--border-light);
      background: #f8f8f8;
      flex-wrap: wrap;
    }

    .batch-toolbar input,
    .batch-toolbar select,
    .params-grid input,
    .manual-range-row input {
      height: 24px;
      padding: 0 6px;
      border: 1px solid var(--border);
      border-radius: 2px;
      background: #fff;
      color: var(--ink);
      font-size: 11px;
    }

    .batch-toolbar input {
      flex: 1;
      min-width: 80px;
    }

    .review-list {
      flex: 1;
      overflow-y: auto;
    }

    .review-row {
      display: flex;
      flex-direction: column;
      gap: 5px;
      width: 100%;
      padding: 7px 8px 8px;
      border: none;
      border-bottom: 1px solid var(--border-inner);
      border-left: 3px solid transparent;
      background: var(--panel);
      text-align: left;
      cursor: pointer;
      min-height: 72px;
    }

    .review-row:hover {
      background: #f0f4f7;
    }

    .review-row[data-range-status="Good"] { border-left-color: var(--good); }
    .review-row[data-range-status="Bad"] { border-left-color: var(--bad); }
    .review-row[data-range-status="Mismatch"],
    .review-row[data-range-status="Error"] { border-left-color: var(--bad); }
    .review-row[data-range-status="Shorter"],
    .review-row[data-range-status="UniprotBetter"] { border-left-color: var(--amber); }

    .review-row-active {
      border-left-color: var(--accent);
      background: #e8eef4;
    }

    .review-row-top {
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 6px;
    }

    .review-row-top strong {
      display: block;
      font-size: 12px;
      line-height: 1.25;
      margin-bottom: 4px;
    }

    .review-row-id {
      font-size: 10px;
      color: var(--muted);
      line-height: 1.35;
      display: -webkit-box;
      -webkit-line-clamp: 2;
      -webkit-box-orient: vertical;
      overflow: hidden;
    }

    .review-row-title {
      min-width: 0;
      flex: 1;
    }

    .review-row-meta {
      display: flex;
      gap: 4px;
      flex-wrap: wrap;
    }

    .review-row-summary {
      display: flex;
      gap: 4px;
      flex-wrap: wrap;
      align-items: center;
    }

    .review-mini-pill {
      display: inline-flex;
      align-items: center;
      min-height: 17px;
      padding: 0 5px;
      border-radius: 999px;
      border: 1px solid var(--border-light);
      background: #f5f5f5;
      font-size: 9px;
      color: #555;
      white-space: nowrap;
    }

    .review-row-dots {
      display: flex;
      gap: 3px;
      align-items: center;
    }

    .status-dot {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-width: 18px;
      height: 15px;
      padding: 0 2px;
      border-radius: 3px;
      font-size: 7px;
      font-weight: 900;
      color: rgba(255, 255, 255, 0.92);
      letter-spacing: 0;
      flex-shrink: 0;
    }

    .dot-green { background: var(--good); }
    .dot-red { background: var(--bad); }
    .dot-amber { background: var(--amber); }
    .dot-gray { background: #aaa; }

    .batch-legend {
      display: flex;
      gap: 8px 12px;
      align-items: center;
      flex-wrap: wrap;
      padding: 5px 8px;
      border-bottom: 1px solid var(--border-inner);
      background: #f4f4f4;
      font-size: 9px;
      color: var(--muted);
    }

    .batch-legend-item {
      display: flex;
      align-items: center;
      gap: 3px;
    }

    .batch-legend-note {
      font-weight: 700;
      color: #555;
    }

    .batch-legend .status-dot {
      min-width: 10px;
      width: 10px;
      height: 10px;
      padding: 0;
      border-radius: 50%;
      font-size: 0;
    }

    .detail-column {
      padding: 10px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      min-width: 0;
    }

    .detail-header {
      display: flex;
      flex-direction: column;
      gap: 6px;
      padding-bottom: 2px;
      border-bottom: 1px solid var(--border-light);
    }

    .detail-header-top {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      align-items: center;
    }

    .detail-header-name {
      font-size: 14px;
      font-weight: 700;
    }

    .detail-evidence-row {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      align-items: center;
    }

    .detail-evidence-pill {
      display: inline-flex;
      align-items: center;
      gap: 5px;
      min-height: 20px;
      padding: 0 7px;
      border-radius: 999px;
      border: 1px solid var(--border-light);
      background: #f5f5f5;
      font-size: 10px;
      color: #333;
      white-space: nowrap;
    }

    .detail-evidence-pill-label {
      font-size: 9px;
      font-weight: 800;
      text-transform: uppercase;
      letter-spacing: 0.03em;
      color: #555;
    }

    .detail-evidence-pill-value {
      font-weight: 600;
      color: #222;
    }

    .pill-gray {
      background: #f0f0f0;
      color: #666;
      border-color: #dbdbdb;
    }

    .metric-chip,
    .pill {
      display: inline-flex;
      align-items: center;
      min-height: 18px;
      padding: 0 6px;
      border-radius: 2px;
      border: 1px solid var(--border-light);
      background: #f5f5f5;
      font-size: 10px;
      white-space: nowrap;
    }

    .pill-green {
      background: rgba(39, 174, 96, 0.12);
      color: #1a7a42;
      border-color: rgba(39, 174, 96, 0.3);
    }

    .pill-red {
      background: rgba(192, 57, 43, 0.1);
      color: #922b21;
      border-color: rgba(192, 57, 43, 0.25);
    }

    .pill-blue {
      background: rgba(40, 104, 160, 0.12);
      color: #20537f;
      border-color: rgba(40, 104, 160, 0.25);
    }

    .pill-amber {
      background: rgba(212, 148, 58, 0.14);
      color: #7a4e10;
      border-color: rgba(212, 148, 58, 0.35);
    }

    .status-ok {
      color: #1a7a42;
      font-weight: 700;
    }

    .status-bad {
      color: #922b21;
      font-weight: 700;
    }

    .status-na {
      color: var(--muted);
      font-weight: 700;
    }

    .detail-section {
      border: 1px solid var(--border-light);
      border-radius: 3px;
      overflow: hidden;
      background: #fcfcfc;
    }

    .metadata-panel {
      display: flex;
      flex-direction: column;
      min-height: calc(100vh - 128px);
      overflow: hidden;
    }

    .metadata-panel-header {
      padding: 9px 10px 7px;
      border-bottom: 1px solid var(--border-light);
      background: #f4f4f4;
    }

    .metadata-panel-header h2 {
      margin: 0;
      font-size: 13px;
    }

    .metadata-panel-header p {
      margin: 3px 0 0;
      font-size: 10px;
      color: var(--muted);
    }

    .metadata-toolbar {
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      align-items: center;
      padding: 8px 10px;
      border-bottom: 1px solid var(--border-light);
      background: #fafafa;
    }

    .metadata-toolbar input {
      flex: 1 1 280px;
      min-width: 220px;
      height: 28px;
      padding: 0 8px;
      border: 1px solid var(--border);
      border-radius: 3px;
      font-size: 12px;
    }

    .metadata-toolbar-note {
      font-size: 11px;
      color: var(--muted);
    }

    .metadata-table-wrap {
      flex: 1;
      overflow: auto;
      background: #fff;
    }

    .metadata-table {
      width: 100%;
      border-collapse: collapse;
      table-layout: auto;
      min-width: 960px;
      font-size: 11px;
    }

    .metadata-table th,
    .metadata-table td {
      padding: 7px 8px;
      border-bottom: 1px solid var(--border-inner);
      border-right: 1px solid #efefef;
      vertical-align: middle;
      text-align: left;
      background: #fff;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }

    .metadata-table th:last-child,
    .metadata-table td:last-child {
      border-right: none;
    }

    .metadata-table thead th {
      position: sticky;
      top: 0;
      z-index: 1;
      padding: 0;
      background: #f1f3f5;
      word-break: normal;
      white-space: nowrap;
    }

    .metadata-th-btn {
      appearance: none;
      width: 100%;
      border: none;
      background: transparent;
      color: #333;
      font: inherit;
      font-size: 11px;
      font-weight: 700;
      text-align: left;
      padding: 8px;
      cursor: pointer;
      white-space: nowrap;
    }

    .metadata-th-btn:hover {
      background: #e8edf1;
    }

    .metadata-table tbody tr:nth-child(even) td {
      background: #fcfcfc;
    }

    .metadata-id-link {
      appearance: none;
      border: none;
      background: transparent;
      color: var(--info);
      font: inherit;
      font-weight: 700;
      padding: 0;
      cursor: pointer;
      text-decoration: underline;
      text-underline-offset: 2px;
      display: block;
      width: 100%;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
    }

    .metadata-id-link:hover {
      color: #1b4f7a;
    }

    .ranges-summary-bar {
      display: flex;
      gap: 6px;
      flex-wrap: wrap;
      align-items: center;
      padding: 8px 10px 0;
      background: #fff;
    }

    .detail-section-header {
      display: flex;
      justify-content: space-between;
      gap: 8px;
      align-items: center;
      padding: 8px 10px;
      background: #f4f4f4;
      border-bottom: 1px solid var(--border-light);
    }

    .detail-section-summary {
      cursor: pointer;
      list-style: none;
    }

    .detail-section-summary::-webkit-details-marker {
      display: none;
    }

    .detail-section-index {
      font-size: 11px;
      font-weight: 700;
    }

    .detail-section-copy {
      font-size: 11px;
      color: var(--muted);
    }

    .section-chip-row {
      display: flex;
      gap: 4px;
      flex-wrap: wrap;
      justify-content: flex-end;
    }

    .track-viewer-wrap {
      padding: 10px;
      background: #fff;
    }

    .detail-top-layout {
      display: grid;
      grid-template-columns: minmax(0, 1.35fr) minmax(320px, 0.9fr);
      gap: 10px;
      align-items: stretch;
    }

    .detail-top-column {
      display: flex;
      flex-direction: column;
      gap: 10px;
      min-width: 0;
    }

    .detail-top-column-right {
      display: flex;
      flex-direction: column;
      min-width: 0;
      min-height: 100%;
      align-self: stretch;
    }

    .detail-top-column-right .detail-section {
      height: 100%;
      display: flex;
      flex-direction: column;
    }

    .detail-top-column-right .metadata-body {
      flex: 1 1 auto;
      min-height: 0;
    }

    .detail-top-column-right .structure-pane-wrap {
      flex: 1 1 auto;
      min-height: 0;
    }

    .detail-top-layout .detail-section {
      margin: 0;
    }

    .structure-docked-host {
      display: flex;
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
      background: #fff;
      border: 1px solid var(--border-light);
      border-radius: 3px;
      overflow: hidden;
    }

    .track-svg {
      display: block;
      width: 100%;
      border: 1px solid var(--border-light);
      background: #fff;
      cursor: crosshair;
      user-select: none;
    }

    .track-svg.dragging {
      cursor: ew-resize;
    }

    .track-drag-note {
      font-size: 11px;
      color: var(--muted);
      margin-top: 6px;
    }

    .coordinates-extra-toggle {
      appearance: none;
      border: 1px solid var(--border);
      background: #fff;
      color: inherit;
      font: inherit;
      font-size: 11px;
      line-height: 1;
      padding: 4px 8px;
      border-radius: 999px;
      cursor: pointer;
    }

    .coordinates-extra-toggle:hover {
      background: #f1f1f1;
    }

    .track-explainer {
      margin-top: 6px;
      border: 1px solid var(--border-light);
      border-radius: 2px;
      background: #fafafa;
      overflow: hidden;
    }

    .track-explainer summary {
      cursor: pointer;
      list-style: none;
      padding: 6px 8px;
      font-size: 11px;
      font-weight: 700;
      background: #f3f3f3;
      border-bottom: 1px solid transparent;
    }

    .track-explainer[open] summary {
      border-bottom-color: var(--border-light);
    }

    .track-explainer summary::-webkit-details-marker {
      display: none;
    }

    .track-explainer-body {
      padding: 8px;
      font-size: 11px;
      color: #444;
      display: flex;
      flex-direction: column;
      gap: 6px;
    }

    .track-explainer-body p {
      margin: 0;
    }

    .track-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 6px;
      font-size: 11px;
      color: var(--muted);
    }

    .track-legend-item {
      display: inline-flex;
      align-items: center;
      gap: 5px;
    }

    .track-legend-swatch {
      display: inline-block;
      width: 18px;
      height: 10px;
      border-radius: 2px;
      border: 1px solid transparent;
      flex: 0 0 auto;
    }

    .track-legend-r1 {
      background: rgba(231, 76, 60, 0.10);
      border-color: #e74c3c;
    }

    .track-legend-r2 {
      background: rgba(230, 126, 34, 0.35);
      border-color: #e67e22;
    }

    .track-legend-r3 {
      background: rgba(39, 174, 96, 0.35);
      border-color: #27ae60;
    }

    .track-legend-custom {
      background: transparent;
      border: 1.5px dashed #7a7a7a;
    }

    .controls-row {
      display: grid;
      grid-template-columns: minmax(240px, 1.4fr) minmax(240px, 1fr);
      gap: 10px;
      padding: 0 10px 10px;
    }

    .controls-section-label {
      font-size: 11px;
      font-weight: 700;
      margin-bottom: 5px;
    }

    .candidate-buttons {
      display: flex;
      gap: 6px;
      flex-wrap: wrap;
    }

    .candidate-btn {
      display: inline-flex;
      align-items: center;
      gap: 6px;
      min-height: 24px;
      padding: 0 8px;
      border: 1px solid var(--border);
      border-radius: 2px;
      background: #fff;
      cursor: pointer;
    }

    .candidate-btn:hover {
      background: #f5f9fc;
    }

    .candidate-btn-active {
      background: #e8eef4;
      border-color: var(--accent);
    }

    .candidate-btn-dot {
      width: 9px;
      height: 9px;
      border-radius: 50%;
      flex: 0 0 auto;
    }

    .manual-range-row,
    .params-grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 8px;
    }

    .manual-range-row {
      margin-bottom: 8px;
    }

    .manual-range-row label,
    .params-field label {
      display: flex;
      flex-direction: column;
      gap: 4px;
      font-size: 11px;
    }

    .range-status-line {
      font-size: 11px;
      color: var(--muted);
    }

    .seq-panel {
      margin-top: 8px;
      display: flex;
      flex-direction: column;
      gap: 4px;
    }

    .seq-track-label {
      font-size: 11px;
      color: var(--muted);
    }

    .seq-pre,
    .output-pre,
    .fasta-pre {
      margin: 0;
      padding: 8px;
      border: 1px solid var(--border-light);
      background: #fff;
      border-radius: 2px;
      overflow-x: auto;
      font-family: var(--font-mono);
      font-size: 12px;
      line-height: 1.5;
    }

    .seq-pre {
      white-space: pre-wrap;
      word-break: break-all;
    }

    .output-pre {
      white-space: pre;
    }

    .fasta-pre {
      white-space: pre-wrap;
      word-break: break-all;
    }

    .seq-ctx {
      color: #777;
    }

    .seq-mark {
      background: #fff0a8;
      color: #000;
      padding: 0;
    }

    .metadata-body {
      padding: 10px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      background: #fff;
    }

    .construct-summary-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
      gap: 8px;
    }

    .summary-card {
      display: flex;
      flex-direction: column;
      gap: 3px;
      min-height: 70px;
      padding: 8px;
      border: 1px solid var(--border-light);
      border-radius: 2px;
      background: #fdfdfd;
    }

    .summary-card span {
      font-size: 10px;
      color: var(--muted);
    }

    .summary-card strong {
      font-size: 12px;
    }

    .metadata-note {
      font-size: 11px;
      color: #444;
      padding: 8px;
      border: 1px dashed var(--border);
      background: #f7f7f7;
    }

    .issue-panel-body {
      padding: 10px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      background: #fff;
    }

    .issue-alert {
      padding: 8px 10px;
      border: 1px solid var(--border-light);
      border-left-width: 4px;
      border-radius: 3px;
      background: #fafafa;
      font-size: 11px;
      color: #333;
    }

    .issue-alert-danger {
      border-left-color: var(--bad);
      background: rgba(192, 57, 43, 0.06);
    }

    .issue-alert-warning {
      border-left-color: var(--amber);
      background: rgba(212, 148, 58, 0.08);
    }

    .issue-subsection {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }

    .issue-subsection-note {
      font-size: 11px;
      color: #555;
    }

    .issue-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px 12px;
      font-size: 10px;
      color: var(--muted);
    }

    .issue-legend-item {
      display: inline-flex;
      align-items: center;
      gap: 5px;
    }

    .issue-legend-swatch {
      width: 14px;
      height: 11px;
      border-radius: 2px;
      border: 1px solid transparent;
      display: inline-block;
      flex: 0 0 auto;
    }

    .issue-legend-swatch-mismatch {
      background: rgba(192, 57, 43, 0.14);
      border-color: rgba(192, 57, 43, 0.4);
    }

    .issue-legend-swatch-missing {
      background: rgba(212, 148, 58, 0.14);
      border-color: rgba(212, 148, 58, 0.4);
    }

    .issue-compare {
      border: 1px solid var(--border-light);
      border-radius: 3px;
      background: #fcfcfc;
      overflow: hidden;
    }

    .issue-compare-block {
      padding: 8px;
      border-top: 1px solid var(--border-inner);
    }

    .issue-compare-block:first-child {
      border-top: none;
    }

    .issue-compare-row {
      display: grid;
      grid-template-columns: 74px 38px minmax(0, 1fr);
      gap: 8px;
      align-items: start;
      font-family: var(--font-mono);
      font-size: 11px;
      line-height: 1.45;
      margin-bottom: 3px;
    }

    .issue-compare-row:last-child {
      margin-bottom: 0;
    }

    .issue-row-label {
      color: #666;
      font-weight: 700;
      text-transform: lowercase;
    }

    .issue-row-start {
      color: #888;
      text-align: right;
    }

    .issue-row-seq,
    .issue-row-cds {
      min-width: 0;
      white-space: normal;
      word-break: break-word;
    }

    .issue-row-seq-diff {
      color: #999;
    }

    .issue-char {
      display: inline-block;
      min-width: 0.72em;
      text-align: center;
      border-radius: 2px;
    }

    .issue-char-mismatch {
      background: rgba(192, 57, 43, 0.14);
      color: #8f2317;
    }

    .issue-char-missing {
      background: rgba(212, 148, 58, 0.14);
      color: #8a5b1f;
    }

    .issue-diff-char {
      display: inline-block;
      min-width: 0.72em;
      text-align: center;
    }

    .issue-diff-mismatch {
      color: #b83a2a;
      font-weight: 700;
    }

    .issue-diff-match {
      color: #b1b1b1;
    }

    .issue-codon {
      display: inline-block;
      min-width: 2.6em;
      margin-right: 0.38em;
      text-align: center;
      border-radius: 2px;
    }

    .issue-codon-mismatch {
      background: rgba(192, 57, 43, 0.14);
      color: #8f2317;
    }

    .issue-inline-missing {
      color: #8a5b1f;
      background: rgba(212, 148, 58, 0.12);
      padding: 0 3px;
      border-radius: 2px;
    }

    .structure-pane-wrap {
      display: flex;
      flex-direction: column;
      align-items: stretch;
      gap: 8px;
      width: 100%;
    }

    .structure-pane-wrap-docked {
      flex: 1 1 auto;
      min-height: 0;
      height: 100%;
    }

    .structure-panel {
      position: relative;
      display: flex;
      flex-direction: column;
      width: 100%;
      height: 460px;
      border: 1px solid var(--border-light);
      background: #fff;
      box-shadow: var(--shadow);
      overflow: hidden;
      flex: 0 0 auto;
    }

    .structure-panel-docked {
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
    }

    .structure-pane-meta {
      display: flex;
      gap: 6px;
      flex-wrap: wrap;
      justify-content: center;
    }

    .structure-viewer {
      position: relative;
      flex: 1 1 auto;
      width: 100%;
      height: 100%;
      min-height: 0;
      background: #fff;
      overflow: hidden;
    }

    .structure-pane-wrap-docked .structure-viewer {
      min-height: 100%;
    }

    .structure-viewer canvas {
      display: block;
      width: 100% !important;
      height: 100% !important;
    }

    .structure-fallback {
      display: flex;
      align-items: center;
      justify-content: center;
      flex: 1 1 auto;
      min-height: 0;
      padding: 16px;
      color: var(--muted);
      font-size: 11px;
      text-align: center;
      border-top: 1px solid var(--border-light);
      background: #fcfcfc;
    }

    .browser-wrap {
      display: flex;
      flex-direction: column;
      gap: 8px;
    }

    .browser-svg {
      display: block;
      width: 100%;
      border: 1px solid var(--border-light);
      background: #fff;
    }

    .browser-empty {
      padding: 12px;
      border: 1px dashed var(--border);
      background: #fafafa;
      color: var(--muted);
      font-size: 11px;
    }

    .browser-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      font-size: 11px;
      color: var(--muted);
    }

    .browser-legend-item {
      display: inline-flex;
      align-items: center;
      gap: 5px;
    }

    .browser-legend-swatch {
      display: inline-block;
      width: 14px;
      height: 10px;
      border-radius: 2px;
      border: 1px solid rgba(0, 0, 0, 0.12);
      flex: 0 0 auto;
    }

    .browser-legend-helix { background: #d94f43; }
    .browser-legend-strand { background: #3a85c8; }
    .browser-legend-coil { background: #c4c4c4; }
    .browser-legend-numeric { background: #6c8ebf; }
    .browser-legend-range { background: rgba(44, 62, 80, 0.14); }

    .domain-chip-row {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
    }

    .domain-chip {
      display: inline-flex;
      align-items: center;
      min-height: 22px;
      padding: 0 8px;
      border: 1px solid var(--border-light);
      background: #f3f7fb;
      border-radius: 2px;
      font-size: 11px;
    }

    .params-panel {
      border-top: 1px solid var(--border-light);
      padding-top: 10px;
    }

    .construct-summary-body {
      padding: 10px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      background: #fff;
    }

    .output-stack {
      display: flex;
      flex-direction: column;
      gap: 8px;
    }

    .empty-state {
      padding: 24px;
      color: var(--muted);
    }

    @media (max-width: 980px) {
      .workspace {
        grid-template-columns: 1fr;
      }

      .detail-top-layout {
        grid-template-columns: 1fr;
      }

      .controls-row {
        grid-template-columns: 1fr;
      }

      .app-toolbar {
        flex-direction: column;
        align-items: stretch;
      }

      .structure-panel {
        height: 360px;
      }
    }
  </style>
</head>
<body>
  <main class="app-shell">
    <div id="toolbar" class="app-toolbar"></div>
    <div id="app-tabs" class="app-tabs"></div>
    <div id="review-page" class="workspace">
      <aside id="batch-panel" class="panel batch-panel"></aside>
      <section id="detail-panel" class="panel detail-column"></section>
    </div>
    <section id="ranges-page" class="panel metadata-panel" hidden></section>
    <section id="metadata-page" class="panel metadata-panel" hidden></section>
  </main>

  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <script id="report-data" type="application/json">__PAYLOAD__</script>
  <script>
    const payload = JSON.parse(document.getElementById("report-data").textContent);
    const dataset = payload.dataset;
    const datasetSummary = payload.datasetSummary;
    const defaultParams = payload.defaultParams;
    const codonTable = payload.codonTable;
    const customRangeIndex = payload.customRanges ?? {};
    const metadataTable = payload.metadataTable ?? null;
    const hasMetadata = Boolean(metadataTable?.columns?.length && metadataTable?.rows?.length);
    const trackPalette = ["#e8b84b", "#3a85c8", "#27ae60", "#9b59b6", "#e74c3c", "#1abc9c", "#d4943a", "#2868a0"];
    const candidateColors = { r1: "#e74c3c", r2: "#e67e22", r3: "#27ae60" };
    const selectedRangeColor = "#111111";
    const ssColors = {
      H: "#d94f43",
      G: "#c43c31",
      I: "#e07a72",
      E: "#3a85c8",
      B: "#2868a0",
      T: "#e8b84b",
      S: "#d4943a",
      P: "#7a8f60",
      "-": "#c4c4c4"
    };
    const ssNames = {
      H: "alpha helix",
      G: "3-10 helix",
      I: "pi helix",
      E: "extended strand",
      B: "isolated beta bridge",
      T: "turn",
      S: "bend",
      P: "polyproline helix",
      "-": "coil / unassigned"
    };
    const structuredCodes = new Set(["H", "G", "I", "E", "B", "T"]);
    const initialEntry =
      dataset.find((entry) => entry.evidence?.structureModel?.text) ??
      dataset.find((entry) => Object.keys(entry.evidence || {}).length > 0) ??
      dataset[0] ??
      null;

    const state = {
      view: "review",
      params: { ...defaultParams },
      search: "",
      filter: "all",
      sort: "id",
      selectedId: initialEntry?.id ?? "",
      manualRanges: {},
      structureDockRight: false,
      rangesSortKey: "id",
      rangesSortDir: "asc",
      metadataSearch: "",
      metadataSortKey: metadataTable?.idKey ?? metadataTable?.columns?.[0]?.key ?? "",
      metadataSortDir: "asc",
    };

    const toolbarEl = document.getElementById("toolbar");
    const tabsEl = document.getElementById("app-tabs");
    const reviewPageEl = document.getElementById("review-page");
    const batchPanelEl = document.getElementById("batch-panel");
    const detailPanelEl = document.getElementById("detail-panel");
    const rangesPageEl = document.getElementById("ranges-page");
    const metadataPageEl = document.getElementById("metadata-page");

    let _structureCache = { entryId: null, viewer: null, containerEl: null, layoutMode: null };
    let _structureInitToken = 0;
    let _analysesCache = { paramsKey: null, result: null };
    let _lastBatchKey = null;
    let _renderTimer = null;

    function captureFocusedTextInput() {
      const active = document.activeElement;
      if (!active || !("id" in active) || !active.id) {
        return null;
      }
      if (!["INPUT", "TEXTAREA"].includes(active.tagName)) {
        return null;
      }
      return {
        id: active.id,
        selectionStart:
          typeof active.selectionStart === "number" ? active.selectionStart : null,
        selectionEnd:
          typeof active.selectionEnd === "number" ? active.selectionEnd : null,
      };
    }

    function restoreFocusedTextInput(snapshot) {
      if (!snapshot?.id) {
        return;
      }
      const next = document.getElementById(snapshot.id);
      if (!next || !["INPUT", "TEXTAREA"].includes(next.tagName)) {
        return;
      }
      next.focus({ preventScroll: true });
      if (
        typeof snapshot.selectionStart === "number" &&
        typeof snapshot.selectionEnd === "number" &&
        typeof next.setSelectionRange === "function"
      ) {
        next.setSelectionRange(snapshot.selectionStart, snapshot.selectionEnd);
      }
    }

    function setPageVisibility() {
      const showReview = state.view === "review";
      const showRanges = state.view === "ranges";
      const showMetadata = state.view === "metadata";
      reviewPageEl.hidden = !showReview;
      rangesPageEl.hidden = !showRanges;
      metadataPageEl.hidden = !showMetadata;
      reviewPageEl.classList.toggle("app-page-hidden", !showReview);
      rangesPageEl.classList.toggle("app-page-hidden", !showRanges);
      metadataPageEl.classList.toggle("app-page-hidden", !showMetadata);
    }

    function debounceRender() {
      clearTimeout(_renderTimer);
      _renderTimer = setTimeout(render, 80);
    }

    function clamp(value, min, max) {
      return Math.max(min, Math.min(max, value));
    }

    function clampRange(range, max) {
      const start = clamp(Number(range.start) || 1, 1, max);
      const end = clamp(Number(range.end) || start, start, max);
      return { start, end };
    }

    function escapeHtml(value) {
      return String(value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#39;");
    }

    function clipDisplayText(value, maxChars = 50) {
      const text = String(value ?? "");
      if (text.length <= maxChars) {
        return text;
      }
      return `${text.slice(0, Math.max(0, maxChars - 3))}...`;
    }

    function metadataSortValue(value) {
      const text = String(value ?? "").trim();
      if (!text) {
        return { kind: "empty", value: "" };
      }
      if (/^[+-]?(?:\d+\.?\d*|\.\d+)$/.test(text)) {
        return { kind: "number", value: Number(text) };
      }
      return { kind: "text", value: text.toLowerCase() };
    }

    function compareMetadataRows(rowA, rowB, key, direction) {
      const left = metadataSortValue(rowA[key] ?? "");
      const right = metadataSortValue(rowB[key] ?? "");
      let order = 0;
      if (left.kind === "empty" && right.kind !== "empty") {
        order = 1;
      } else if (left.kind !== "empty" && right.kind === "empty") {
        order = -1;
      } else if (left.kind === "number" && right.kind === "number") {
        order = left.value - right.value;
      } else {
        order = String(left.value).localeCompare(String(right.value), undefined, { numeric: true });
      }
      if (order === 0 && metadataTable?.idKey && key !== metadataTable.idKey) {
        order = String(rowA[metadataTable.idKey] ?? "").localeCompare(String(rowB[metadataTable.idKey] ?? ""));
      }
      return direction === "desc" ? -order : order;
    }

    function compareRangeTableRows(rowA, rowB, key, direction) {
      function compareNumbers(left, right) {
        if (left === "" && right !== "") return 1;
        if (left !== "" && right === "") return -1;
        if (left === "" && right === "") return 0;
        return Number(left) - Number(right);
      }

      let order = 0;
      if (key === "aaRange") {
        order = compareNumbers(rowA.aaStart, rowB.aaStart) || compareNumbers(rowA.aaEnd, rowB.aaEnd);
      } else if (key === "customRange") {
        order = compareNumbers(rowA.customStart, rowB.customStart) || compareNumbers(rowA.customEnd, rowB.customEnd);
      } else if (key === "automaticRange") {
        order = compareNumbers(rowA.automaticStart, rowB.automaticStart) || compareNumbers(rowA.automaticEnd, rowB.automaticEnd);
      } else if (key === "cdsRange") {
        order = compareNumbers(rowA.cdsStart, rowB.cdsStart) || compareNumbers(rowA.cdsEnd, rowB.cdsEnd);
      } else if (["aaLength", "cdsLength", "aaStart", "aaEnd", "cdsStart", "cdsEnd"].includes(key)) {
        order = compareNumbers(rowA[key] ?? "", rowB[key] ?? "");
      } else {
        order = compareMetadataRows(rowA, rowB, key, "asc");
      }

      if (order === 0 && key !== "id") {
        order = String(rowA.id ?? "").localeCompare(String(rowB.id ?? ""), undefined, { numeric: true });
      }
      return direction === "desc" ? -order : order;
    }

    tabsEl.addEventListener("click", (event) => {
      const button = event.target.closest("[data-app-view]");
      if (!button) {
        return;
      }
      const nextView = button.getAttribute("data-app-view") || "review";
      if (state.view === nextView) {
        return;
      }
      state.view = nextView;
      render();
    });

    function renderTabs() {
      tabsEl.style.display = "flex";
      tabsEl.innerHTML = `
        <button type="button" class="app-tab ${state.view === "review" ? "app-tab-active" : ""}" data-app-view="review">
          Review
        </button>
        <button type="button" class="app-tab ${state.view === "ranges" ? "app-tab-active" : ""}" data-app-view="ranges">
          Ranges (${dataset.length})
        </button>
        ${
          hasMetadata
            ? `<button type="button" class="app-tab ${state.view === "metadata" ? "app-tab-active" : ""}" data-app-view="metadata">
                 Metadata (${metadataTable.rows.length})
               </button>`
            : ""
        }
      `;
      if (state.view === "metadata" && !hasMetadata) {
        state.view = "review";
      }
      setPageVisibility();
    }

    function renderRangesTab(analyses) {
      const rows = buildRangeSelectionRows(analyses);
      const sortKey = state.rangesSortKey || "id";
      const sortedRows = rows.slice().sort((a, b) => compareRangeTableRows(a, b, sortKey, state.rangesSortDir));
      const exportHref = buildExportHref(buildRangesExportTsv(sortedRows));
      const okCount = rows.filter((row) => row.constructStatus === "OK").length;
      const manualCount = rows.filter((row) => row.source === "manual").length;
      const customDiffCount = rows.filter((row) => row.customAuto === "different").length;
      const rangeColumns = [
        { key: "id", label: "ID", width: "28ch" },
        { key: "source", label: "source", width: "10ch" },
        { key: "aaRange", label: "selected range", width: "12ch" },
        { key: "customRange", label: "custom range", width: "13ch" },
        { key: "automaticRange", label: "automatic range", width: "13ch" },
        { key: "customAuto", label: "Custom vs Automatic", width: "18ch" },
        { key: "aaLength", label: "AA len", width: "12ch" },
        { key: "cdsRange", label: "CDS range", width: "12ch" },
        { key: "cdsLength", label: "CDS len", width: "11ch" },
        { key: "constructStatus", label: "construct", width: "12ch" },
        { key: "automaticLabel", label: "auto ref", width: "13ch" },
      ];

      rangesPageEl.innerHTML = `
        <div class="metadata-panel-header">
          <h2>Selected Ranges</h2>
          <p>Current per-protein selected ranges. Manual edits in the review pane are reflected here immediately.</p>
        </div>
        <div class="metadata-toolbar">
          <span class="metadata-toolbar-note">${rows.length} proteins · ${okCount} export-ready · ${manualCount} manual · ${customDiffCount} custom/automatic differences</span>
          <a class="action-button action-button-secondary" href="${exportHref}" download="selected_ranges.tsv">Export ranges TSV</a>
        </div>
        <div class="ranges-summary-bar">
          <span class="metric-chip">rows: ${rows.length}</span>
          <span class="metric-chip">export-ready: ${okCount}</span>
          <span class="metric-chip">manual: ${manualCount}</span>
          <span class="metric-chip ${customDiffCount ? "pill-amber" : "pill-green"}">custom != automatic: ${customDiffCount}</span>
          <span class="metric-chip">sorted by: ${escapeHtml(sortKey)} ${state.rangesSortDir}</span>
        </div>
        <div class="metadata-table-wrap">
          <table class="metadata-table">
            <colgroup>
              ${rangeColumns.map((column) => `<col style="width:${column.width}">`).join("")}
            </colgroup>
            <thead>
              <tr>
                ${rangeColumns.map((column) => {
                  const isActive = state.rangesSortKey === column.key;
                  const arrow = isActive ? (state.rangesSortDir === "asc" ? " ↑" : " ↓") : "";
                  return `
                    <th style="min-width:${column.width}">
                      <button type="button" class="metadata-th-btn" data-range-sort="${escapeHtml(column.key)}">
                        ${escapeHtml(column.label)}${arrow}
                      </button>
                    </th>
                  `;
                }).join("")}
              </tr>
            </thead>
            <tbody>
              ${sortedRows.map((row) => `
                <tr>
                  <td style="min-width:28ch" title="${escapeHtml(row.id)}"><button type="button" class="metadata-id-link" data-range-select-id="${escapeHtml(row.id)}" title="${escapeHtml(row.id)}">${escapeHtml(row.id)}</button></td>
                  <td style="min-width:10ch" title="${escapeHtml(row.source)}">${
                    row.source === "manual"
                      ? '<span class="review-mini-pill pill-amber">manual</span>'
                      : escapeHtml(row.source)
                  }</td>
                  <td style="min-width:12ch" title="${escapeHtml(row.aaRange)}">${escapeHtml(row.aaRange)}</td>
                  <td style="min-width:13ch" title="${escapeHtml(row.customRange || 'NA')}">${row.customRange ? escapeHtml(row.customRange) : '<span class="status-na">NA</span>'}</td>
                  <td style="min-width:13ch" title="${escapeHtml(row.automaticRange || 'NA')}">${
                    row.automaticRange
                      ? escapeHtml(row.automaticRange)
                      : '<span class="status-na">NA</span>'
                  }</td>
                  <td style="min-width:10ch">${
                    row.customAuto === "different"
                      ? `<span class="review-mini-pill pill-amber" title="${escapeHtml(`custom ${row.customRange} differs from automatic ${row.automaticRange}`)}">different</span>`
                      : row.customAuto === "match"
                        ? '<span class="review-mini-pill pill-green">match</span>'
                        : '<span class="status-na">NA</span>'
                  }</td>
                  <td style="min-width:12ch" title="${escapeHtml(String(row.aaLength))}">${escapeHtml(String(row.aaLength))}</td>
                  <td style="min-width:12ch" title="${escapeHtml(row.cdsRange || "NA")}">${escapeHtml(row.cdsRange || "NA")}</td>
                  <td style="min-width:11ch" title="${escapeHtml(String(row.cdsLength || "NA"))}">${escapeHtml(String(row.cdsLength || "NA"))}</td>
                  <td style="min-width:12ch"><span class="review-mini-pill ${row.constructStatus === "OK" ? "pill-green" : "pill-red"}">${escapeHtml(row.constructStatus)}</span></td>
                  <td style="min-width:13ch" title="${escapeHtml(row.automaticLabel)}">${escapeHtml(row.automaticLabel)}</td>
                </tr>
              `).join("")}
            </tbody>
          </table>
        </div>
      `;

      rangesPageEl.querySelectorAll("[data-range-sort]").forEach((button) => {
        button.addEventListener("click", () => {
          const key = button.getAttribute("data-range-sort") || "";
          if (!key) {
            return;
          }
          if (state.rangesSortKey === key) {
            state.rangesSortDir = state.rangesSortDir === "asc" ? "desc" : "asc";
          } else {
            state.rangesSortKey = key;
            state.rangesSortDir = "asc";
          }
          render();
        });
      });

      rangesPageEl.querySelectorAll("[data-range-select-id]").forEach((button) => {
        button.addEventListener("click", () => {
          const targetId = button.getAttribute("data-range-select-id");
          if (!targetId) {
            return;
          }
          state.selectedId = targetId;
          state.view = "review";
          render();
        });
      });
    }

    function renderMetadataTab() {
      if (!hasMetadata) {
        metadataPageEl.innerHTML = `<div class="empty-state">No metadata table is available for this report.</div>`;
        return;
      }

      const columns = metadataTable.columns;
      const columnSpecs = columns.map((column) => ({
        ...column,
        minWidthCh: Math.max(10, String(column.label || column.key || "").length + 2),
      }));
      const searchNeedle = state.metadataSearch.trim().toLowerCase();
      const filteredRows = metadataTable.rows.filter((row) => {
        if (!searchNeedle) {
          return true;
        }
        return columnSpecs.some((column) => String(row[column.key] ?? "").toLowerCase().includes(searchNeedle));
      });
      const sortKey = state.metadataSortKey || metadataTable.idKey || columnSpecs[0].key;
      const sortedRows = filteredRows
        .slice()
        .sort((a, b) => compareMetadataRows(a, b, sortKey, state.metadataSortDir));

      metadataPageEl.innerHTML = `
        <div class="metadata-panel-header">
          <h2>Metadata Table</h2>
          <p>${escapeHtml(metadataTable.source || "Optional metadata")} · rows mapped to proteins in this report.</p>
        </div>
        <div class="metadata-toolbar">
          <input
            id="metadata-search-input"
            value="${escapeHtml(state.metadataSearch)}"
            placeholder="Search metadata across all columns…"
            aria-label="Search metadata"
          >
          <span class="metadata-toolbar-note">
            showing ${sortedRows.length} / ${metadataTable.rows.length} rows · sorted by ${escapeHtml(sortKey)} ${state.metadataSortDir}
          </span>
        </div>
        <div class="metadata-table-wrap">
          ${
            sortedRows.length
              ? `<table class="metadata-table">
                  <colgroup>
                    ${columnSpecs.map((column) => `<col style="width:${column.minWidthCh}ch">`).join("")}
                  </colgroup>
                  <thead>
                    <tr>
                      ${columnSpecs.map((column) => {
                        const isActive = state.metadataSortKey === column.key;
                        const arrow = isActive ? (state.metadataSortDir === "asc" ? " ↑" : " ↓") : "";
                        return `
                          <th style="min-width:${column.minWidthCh}ch">
                            <button type="button" class="metadata-th-btn" data-metadata-sort="${escapeHtml(column.key)}">
                              ${escapeHtml(column.label)}${arrow}
                            </button>
                          </th>
                        `;
                      }).join("")}
                    </tr>
                  </thead>
                  <tbody>
                    ${sortedRows.map((row) => `
                      <tr>
                        ${columnSpecs.map((column) => {
                          const value = String(row[column.key] ?? "");
                          if (column.key === metadataTable.idKey && dataset.some((entry) => entry.id === value)) {
                            return `<td style="min-width:${column.minWidthCh}ch" title="${escapeHtml(value)}"><button type="button" class="metadata-id-link" data-metadata-select-id="${escapeHtml(value)}" title="${escapeHtml(value)}">${escapeHtml(clipDisplayText(value, 50))}</button></td>`;
                          }
                          return `<td style="min-width:${column.minWidthCh}ch" title="${escapeHtml(value)}">${escapeHtml(clipDisplayText(value, 50))}</td>`;
                        }).join("")}
                      </tr>
                    `).join("")}
                  </tbody>
                </table>`
              : `<div class="empty-state">No metadata rows matched the current search.</div>`
          }
        </div>
      `;

      metadataPageEl.querySelector("#metadata-search-input")?.addEventListener("input", (event) => {
        state.metadataSearch = event.target.value;
        debounceRender();
      });
      metadataPageEl.querySelectorAll("[data-metadata-sort]").forEach((button) => {
        button.addEventListener("click", () => {
          const key = button.getAttribute("data-metadata-sort") || "";
          if (!key) {
            return;
          }
          if (state.metadataSortKey === key) {
            state.metadataSortDir = state.metadataSortDir === "asc" ? "desc" : "asc";
          } else {
            state.metadataSortKey = key;
            state.metadataSortDir = "asc";
          }
          render();
        });
      });
      metadataPageEl.querySelectorAll("[data-metadata-select-id]").forEach((button) => {
        button.addEventListener("click", () => {
          const targetId = button.getAttribute("data-metadata-select-id");
          if (!targetId) {
            return;
          }
          state.selectedId = targetId;
          state.view = "review";
          render();
        });
      });
    }

    function sanitizeDomId(value) {
      return String(value).replace(/[^A-Za-z0-9_-]+/g, "-");
    }

    function formatRange(range) {
      if (!range) {
        return "NA";
      }
      return `${range.start}-${range.end}`;
    }

    function parseRangeString(value) {
      const [start, end] = String(value || "").split("-").map(Number);
      if (!start || !end) {
        return null;
      }
      return { start, end };
    }

    function toCdsRange(range) {
      return {
        start: (range.start - 1) * 3 + 1,
        end: range.end * 3
      };
    }

    function translateDna(sequence) {
      const normalized = String(sequence || "").toUpperCase();
      let output = "";

      for (let index = 0; index + 2 < normalized.length; index += 3) {
        const codon = normalized.slice(index, index + 3);
        output += codonTable[codon] ?? "X";
      }

      return output;
    }

    function validateProtein(entry) {
      const proteinLength = entry.proteinSequence.length;
      const cdsLength = entry.cdsSequence.length;
      let lengthStatus = "Error";

      if (proteinLength * 3 === cdsLength) {
        lengthStatus = "Match";
      } else if (proteinLength * 3 === cdsLength - 3) {
        lengthStatus = "MatchSTOP";
      }

      const translated = translateDna(entry.cdsSequence).replace(/\*$/, "");
      const translationStatus = translated === entry.proteinSequence ? "Match" : "Mismatch";

      // Find first divergence position for diagnostic detail
      let firstMismatchAa = null;
      if (translationStatus === "Mismatch") {
        const minLen = Math.min(translated.length, entry.proteinSequence.length);
        for (let i = 0; i < minLen; i++) {
          if (translated[i] !== entry.proteinSequence[i]) {
            firstMismatchAa = i + 1;
            break;
          }
        }
        if (firstMismatchAa === null) firstMismatchAa = minLen + 1;
      }

      const cdsPeptideLength = translated.length;
      const lengthDelta = cdsPeptideLength - proteinLength;

      return {
        lengthStatus,
        translationStatus,
        firstMismatchAa,
        cdsPeptideLength,
        lengthDelta,
        exportReady:
          (lengthStatus === "Match" || lengthStatus === "MatchSTOP") &&
          translationStatus === "Match"
      };
    }

    function buildConstruct(entry, range) {
      const cdsRange = toCdsRange(range);
      if (cdsRange.end > entry.cdsSequence.length) {
        return {
          range,
          cdsRange: null,
          peptide: "",
          cds: "",
          translation: "",
          status: "Range exceeds CDS length"
        };
      }

      const peptide = entry.proteinSequence.slice(range.start - 1, range.end);
      const cds = entry.cdsSequence.slice(cdsRange.start - 1, cdsRange.end);
      const translation = translateDna(cds).replace(/\*$/, "");

      return {
        range,
        cdsRange,
        peptide,
        cds,
        translation,
        status: translation === peptide ? "OK" : "Translation mismatch"
      };
    }

    function normalizeRange(range, proteinLength) {
      const start = clamp(range.start, 1, proteinLength);
      const end = clamp(range.end, start, proteinLength);
      return { start, end };
    }

    function mergeDomainsAcrossProtein(domains) {
      if (!domains.length) {
        return [];
      }

      const sortedDomains = domains
        .slice()
        .sort((a, b) => a.start - b.start || a.end - b.end || a.label.localeCompare(b.label));
      const labels = [];
      sortedDomains.forEach((domain) => {
        if (domain.label && !labels.includes(domain.label)) {
          labels.push(domain.label);
        }
      });

      return [{
        start: Math.min(...sortedDomains.map((domain) => domain.start)),
        end: Math.max(...sortedDomains.map((domain) => domain.end)),
        label: labels.join(",") || "merged"
      }];
    }

    function effectiveMergedDomains(entry) {
      if (entry.individualDomains.length) {
        return mergeDomainsAcrossProtein(entry.individualDomains);
      }
      if (entry.mergedDomains.length) {
        return mergeDomainsAcrossProtein(entry.mergedDomains);
      }
      return [];
    }

    function envelopeRange(domains) {
      if (!domains.length) {
        return null;
      }

      return {
        start: Math.min(...domains.map((domain) => domain.start)),
        end: Math.max(...domains.map((domain) => domain.end))
      };
    }

    function expandDomainBlocks(domains, slop, proteinLength) {
      if (!domains.length) {
        return [];
      }

      return mergeDomainsAcrossProtein(
        domains.map((domain) => ({
          ...domain,
          start: Math.max(1, domain.start - slop),
          end: Math.min(proteinLength, domain.end + slop)
        })),
      );
    }

    function preferredStructureTrack(entry) {
      if (entry.evidence.structureDssp?.compatible) {
        return entry.evidence.structureDssp;
      }
      if (entry.evidence.structureUniprot?.compatible) {
        return entry.evidence.structureUniprot;
      }
      return null;
    }

    function buildRuns(values) {
      const runs = [];
      let current = null;

      values.forEach((value, index) => {
        if (current && current.value === value) {
          current.end = index;
        } else {
          if (current) {
            runs.push(current);
          }
          current = { value, start: index, end: index };
        }
      });

      if (current) {
        runs.push(current);
      }

      return runs;
    }

    function structureTypeName(code) {
      if (code === null || code === undefined) {
        return "no assignment";
      }
      return ssNames[code] ?? `other (${code})`;
    }

    function getStructureSelection(entry, activeRange) {
      const model = entry.evidence.structureModel;
      const mapping = model?.mapping;
      if (!model?.text || !mapping) {
        return {
          hasMapping: false,
          hasOverlap: false,
          proteinRange: null,
          residueRange: null
        };
      }

      const proteinStart = Number(mapping.proteinStart);
      const proteinEnd = Number(mapping.proteinEnd);
      const residueNumbers = Array.isArray(mapping.residueNumbers) ? mapping.residueNumbers : [];
      if (!proteinStart || !proteinEnd || !residueNumbers.length) {
        return {
          hasMapping: false,
          hasOverlap: false,
          proteinRange: null,
          residueRange: null
        };
      }

      const overlapStart = Math.max(activeRange.start, proteinStart);
      const overlapEnd = Math.min(activeRange.end, proteinEnd);
      if (overlapStart > overlapEnd) {
        return {
          hasMapping: true,
          hasOverlap: false,
          proteinRange: { start: proteinStart, end: proteinEnd },
          residueRange: null,
          mappingSource: mapping.mappingSource ?? null
        };
      }

      const localStartIndex = overlapStart - proteinStart;
      const localEndIndex = overlapEnd - proteinStart;
      const residueStart = residueNumbers[localStartIndex] ?? null;
      const residueEnd = residueNumbers[localEndIndex] ?? null;

      return {
        hasMapping: true,
        hasOverlap: residueStart !== null && residueEnd !== null,
        proteinRange: { start: proteinStart, end: proteinEnd },
        overlapRange: { start: overlapStart, end: overlapEnd },
        residueRange:
          residueStart !== null && residueEnd !== null
            ? { start: residueStart, end: residueEnd }
            : null,
        mappingSource: mapping.mappingSource ?? null
      };
    }

    function buildAreaPath(values, xPos, top, height, maxValue = 1) {
      const bottom = top + height;
      const segments = [];
      let currentPoints = [];

      const flush = () => {
        if (currentPoints.length >= 2) {
          const firstX = currentPoints[0].split(",")[0];
          const lastX = currentPoints[currentPoints.length - 1].split(",")[0];
          segments.push(`M ${firstX} ${bottom} L ${currentPoints.join(" L ")} L ${lastX} ${bottom} Z`);
        }
        currentPoints = [];
      };

      values.forEach((value, index) => {
        if (value === null || value === undefined) {
          flush();
          return;
        }
        const x = xPos(index + 1);
        const normalized = Math.max(0, Math.min(1, Number(value) / maxValue));
        const y = bottom - normalized * height;
        currentPoints.push(`${x.toFixed(1)},${y.toFixed(1)}`);
      });

      flush();
      return segments.join(" ");
    }

    function extractStructuredRuns(track, minRunLength) {
      const runs = [];
      let currentStart = null;

      track.values.forEach((value, index) => {
        const isStructured = value ? structuredCodes.has(value) : false;
        if (isStructured && currentStart === null) {
          currentStart = index + 1;
        }
        if (!isStructured && currentStart !== null) {
          const run = { start: currentStart, end: index };
          if (run.end - run.start + 1 >= minRunLength) {
            runs.push(run);
          }
          currentStart = null;
        }
      });

      if (currentStart !== null) {
        const run = { start: currentStart, end: track.values.length };
        if (run.end - run.start + 1 >= minRunLength) {
          runs.push(run);
        }
      }

      return runs;
    }

    function expandRangeFromNearbyStructure(baseRange, structuredRuns, offset, proteinLength) {
      const next = { ...baseRange };

      structuredRuns.forEach((run) => {
        const touchesLeft = run.start <= next.start + offset && run.end >= next.start - offset;
        const touchesRight = run.start <= next.end + offset && run.end >= next.end - offset;

        if (touchesLeft) {
          next.start = Math.min(next.start, Math.max(1, run.start - 1));
        }
        if (touchesRight) {
          next.end = Math.max(next.end, Math.min(proteinLength, run.end + 1));
        }
      });

      return normalizeRange(next, proteinLength);
    }

    function buildSuggestions(entry, params) {
      const r1Segments = effectiveMergedDomains(entry);
      const domainRange = envelopeRange(r1Segments);
      if (!domainRange) {
        return { candidates: {}, recommendedKey: null };
      }

      const proteinLength = entry.proteinSequence.length;
      const r1 = {
        ...normalizeRange(domainRange, proteinLength),
        segments: r1Segments
      };
      const r2Segments = expandDomainBlocks(r1Segments, params.slop, proteinLength);
      const r2 = {
        ...normalizeRange(envelopeRange(r2Segments) ?? domainRange, proteinLength),
        segments: r2Segments
      };

      const structure = preferredStructureTrack(entry);
      const structuredRuns =
        structure && structure.compatible
          ? extractStructuredRuns(structure, params.minStructuredRun)
          : [];

      let r3 = expandRangeFromNearbyStructure(
        r2,
        structuredRuns,
        params.offset,
        proteinLength
      );

      if (r3.start > 1 && r3.start < params.nTerminalSnapThreshold) {
        r3 = { ...r3, start: 1 };
      }
      if (r3.end < proteinLength && proteinLength - r3.end < params.nTerminalSnapThreshold) {
        r3 = { ...r3, end: proteinLength };
      }

      return {
        candidates: {
          r1: {
            key: "r1",
            label: "r1",
            description: entry.individualDomains.length
              ? "All individual domain hits merged into one protein-level span"
              : "All domain intervals merged into one protein-level span",
            ...r1
          },
          r2: {
            key: "r2",
            label: "r2",
            description: `r1 plus ±${params.slop} aa`,
            ...r2
          },
          r3: {
            key: "r3",
            label: "r3",
            description: "Structure-aware expansion around r2",
            ...r3
          }
        },
        recommendedKey: "r3"
      };
    }

    function buildAnalyses() {
      const paramsKey = JSON.stringify(state.params);
      if (_analysesCache.paramsKey === paramsKey && _analysesCache.result) {
        return _analysesCache.result;
      }
      function normalizeCustomRange(rangeLike, proteinLength) {
        if (!rangeLike) {
          return null;
        }

        let parsed = null;
        if (typeof rangeLike === "string") {
          parsed = parseRangeString(rangeLike);
        } else if (typeof rangeLike === "object") {
          if (Number(rangeLike.start) && Number(rangeLike.end)) {
            parsed = { start: Number(rangeLike.start), end: Number(rangeLike.end) };
          } else if (rangeLike.range) {
            parsed = parseRangeString(rangeLike.range);
          }
        }

        if (!parsed) {
          return null;
        }

        const normalized = clampRange(parsed, proteinLength);
        return {
          start: normalized.start,
          end: normalized.end,
          label:
            (typeof rangeLike === "object" && (rangeLike.label || rangeLike.source)) ||
            "custom range"
        };
      }

      const result = dataset.map((entry) => ({
        entry,
        suggestion: buildSuggestions(entry, state.params),
        validation: validateProtein(entry),
        referenceRange: parseRangeString(entry.reference?.picked_range ?? ""),
        customRanges: (Array.isArray(customRangeIndex[entry.id]) ? customRangeIndex[entry.id] : [])
          .map((rangeLike) => normalizeCustomRange(rangeLike, entry.proteinSequence.length))
          .filter(Boolean),
      }));
      _analysesCache = { paramsKey, result };
      return result;
    }

    function primaryCustomRange(analysis) {
      return analysis.customRanges.length ? analysis.customRanges[0] : null;
    }

    function automaticRange(analysis) {
      return (
        (analysis.suggestion.recommendedKey
          ? analysis.suggestion.candidates[analysis.suggestion.recommendedKey]
          : null) ?? analysis.referenceRange ?? null
      );
    }

    function rangesEqual(left, right) {
      return Boolean(left && right && left.start === right.start && left.end === right.end);
    }

    function customAutoDifference(analysis) {
      const custom = primaryCustomRange(analysis);
      const automatic = automaticRange(analysis);
      if (!custom || !automatic) {
        return null;
      }
      return {
        custom,
        automatic,
        different: !rangesEqual(custom, automatic),
        automaticLabel: analysis.suggestion.recommendedKey ?? analysis.entry.reference?.picked_range_name ?? "auto",
      };
    }

    function analysisHasCdsProteinProblem(analysis) {
      return (
        analysis.validation.translationStatus !== "Match" ||
        !["Match", "MatchSTOP"].includes(analysis.validation.lengthStatus)
      );
    }

    function analysisReviewStatus(analysis) {
      if (analysisHasCdsProteinProblem(analysis)) {
        return "Mismatch";
      }
      return analysis.entry.reference?.status_range ?? "";
    }

    function analysisGroup(analysis) {
      if (analysisHasCdsProteinProblem(analysis)) {
        return "CDS error";
      }
      const customDiff = customAutoDifference(analysis);
      if (customDiff?.different) {
        return "Range mismatch";
      }
      const reviewStatus = analysisReviewStatus(analysis);
      if (reviewStatus && reviewStatus !== "Good") {
        return "Range mismatch";
      }
      return "Good";
    }

    function defaultRangeSourceLabel(analysis) {
      if (primaryCustomRange(analysis)) {
        return "custom";
      }
      return analysis.suggestion.recommendedKey ?? analysis.entry.reference?.picked_range_name ?? "range";
    }

    function getActiveRange(analysis) {
      const customRange = primaryCustomRange(analysis);
      const recommendedRange = automaticRange(analysis);

      return clampRange(
        state.manualRanges[analysis.entry.id] ??
          customRange ??
          recommendedRange ?? {
            start: 1,
            end: Math.min(50, analysis.entry.proteinSequence.length)
          },
        analysis.entry.proteinSequence.length
      );
    }

    function dotClass(value) {
      if (["Good", "Full", "MatchSTOP", "Match"].includes(value)) {
        return "dot-green";
      }
      if (["Bad", "Missing", "Error", "Mismatch"].includes(value)) {
        return "dot-red";
      }
      if (["Shorter", "UniprotBetter"].includes(value)) {
        return "dot-amber";
      }
      return "dot-gray";
    }

    function evidenceBadges(entry) {
      const ev = entry.evidence ?? {};
      const has3d = Boolean(ev.structureModel?.text);
      const hasDssp = ev.structureDssp?.compatible;
      const hasUniprot = ev.structureUniprot?.compatible;
      const hasCons = ev.conservation?.compatible;
      const hasConsFull = ev.conservationFull?.compatible;
      const hasIupred = ev.iupred?.compatible;
      const hasPlddt = ev.plddt?.compatible;

      let ssTitle;
      if (hasDssp) {
        ssTitle = `secondary structure: DSSP (${Math.round((ev.structureDssp.coverage ?? 1) * 100)}% coverage)`;
      } else if (hasUniprot) {
        ssTitle = "secondary structure: UniProt only";
      } else if (has3d) {
        ssTitle = "secondary structure: absent — DSSP not computed for this protein (PDB available)";
      } else {
        ssTitle = "secondary structure: absent";
      }

      return [
        {
          letter: "3D",
          title: has3d ? "3D structure: present" : "3D structure: absent",
          cls: has3d ? "dot-green" : "dot-gray",
          score: has3d ? 1 : 0,
        },
        {
          letter: "SS",
          title: ssTitle,
          cls: hasDssp
            ? ((ev.structureDssp.coverage ?? 1) < 1 ? "dot-amber" : "dot-green")
            : hasUniprot ? "dot-amber" : "dot-gray",
          score: hasDssp ? 2 : hasUniprot ? 1 : 0,
        },
        {
          letter: "CN",
          title: (hasCons && hasConsFull) ? "conservation: DBD + full" : hasConsFull ? "conservation: full only" : hasCons ? "conservation: DBD only" : "conservation: absent",
          cls: (hasCons || hasConsFull) ? "dot-green" : "dot-gray",
          score: (hasCons ? 1 : 0) + (hasConsFull ? 1 : 0),
        },
        {
          letter: "ID",
          title: hasIupred ? "IUPred disorder: present" : "IUPred disorder: absent",
          cls: hasIupred ? "dot-green" : "dot-gray",
          score: hasIupred ? 1 : 0,
        },
        {
          letter: "pL",
          title: hasPlddt ? "pLDDT: present" : "pLDDT: absent",
          cls: hasPlddt ? "dot-green" : "dot-gray",
          score: hasPlddt ? 1 : 0,
        },
      ];
    }

    function evidenceCards(entry, validation) {
      const ev = entry.evidence ?? {};
      const has3d = Boolean(ev.structureModel?.text);
      const hasDssp = ev.structureDssp?.compatible;
      const hasUniprot = ev.structureUniprot?.compatible;
      const hasCons = ev.conservation?.compatible;
      const hasConsFull = ev.conservationFull?.compatible;
      const hasIupred = ev.iupred?.compatible;
      const hasPlddt = ev.plddt?.compatible;
      const dsspCoverage = Math.round((ev.structureDssp?.coverage ?? 1) * 100);

      const cdsCard = (() => {
        if (validation.exportReady) {
          return {
            label: "CDS",
            value: "match",
            cls: "dot-green",
            title: "CDS: translation matches protein"
          };
        }
        if (validation.translationStatus === "Mismatch") {
          const parts = ["CDS: translation mismatch"];
          if (validation.firstMismatchAa !== null) {
            parts.push(`first divergence at aa ${validation.firstMismatchAa}`);
          }
          if (validation.lengthDelta !== 0) {
            parts.push(
              `CDS encodes ${validation.cdsPeptideLength} aa (${validation.lengthDelta > 0 ? "+" : ""}${validation.lengthDelta} vs protein)`
            );
          }
          return {
            label: "CDS",
            value: "mismatch",
            cls: "dot-red",
            title: parts.join(", ")
          };
        }
        return {
          label: "CDS",
          value: validation.lengthStatus.toLowerCase(),
          cls: "dot-amber",
          title: `CDS: length ${validation.lengthStatus}`
        };
      })();

      return [
        {
          label: "3D model",
          value: has3d ? "present" : "absent",
          cls: has3d ? "dot-green" : "dot-gray",
          title: has3d ? "3D structure: present" : "3D structure: absent",
        },
        {
          label: "2° structure",
          value: hasDssp ? `DSSP ${dsspCoverage}%` : hasUniprot ? "UniProt only" : "absent",
          cls: hasDssp ? (dsspCoverage < 100 ? "dot-amber" : "dot-green") : hasUniprot ? "dot-amber" : "dot-gray",
          title: hasDssp
            ? `secondary structure: DSSP (${dsspCoverage}% coverage)`
            : hasUniprot
              ? "secondary structure: UniProt only"
              : has3d
                ? "secondary structure: absent — DSSP not computed for this protein (PDB available)"
                : "secondary structure: absent",
        },
        {
          label: "Conservation",
          value: hasCons && hasConsFull ? "full + DBD" : hasConsFull ? "full only" : hasCons ? "DBD only" : "absent",
          cls: hasCons || hasConsFull ? "dot-green" : "dot-gray",
          title: hasCons && hasConsFull
            ? "conservation: DBD + full"
            : hasConsFull
              ? "conservation: full only"
              : hasCons
                ? "conservation: DBD only"
                : "conservation: absent",
        },
        {
          label: "IUPred",
          value: hasIupred ? "present" : "absent",
          cls: hasIupred ? "dot-green" : "dot-gray",
          title: hasIupred ? "IUPred disorder: present" : "IUPred disorder: absent",
        },
        {
          label: "pLDDT",
          value: hasPlddt ? "present" : "absent",
          cls: hasPlddt ? "dot-green" : "dot-gray",
          title: hasPlddt ? "pLDDT: present" : "pLDDT: absent",
        },
        cdsCard,
      ];
    }

    function evidenceSortScore(entry, key) {
      const badges = evidenceBadges(entry);
      const byLetter = Object.fromEntries(badges.map((b) => [b.letter, b.score]));
      if (key === "structure") return -(byLetter["3D"] * 4 + byLetter["SS"]);
      if (key === "conservation") return -byLetter["CN"];
      if (key === "disorder") return -byLetter["ID"];
      if (key === "plddt") return -byLetter["pL"];
      if (key === "evidence") return -(byLetter["3D"] + byLetter["SS"] + byLetter["CN"] + byLetter["ID"] + byLetter["pL"]);
      return 0;
    }

    function dotToTone(value) {
      const dc = dotClass(value);
      if (dc === "dot-green") return "pill-green";
      if (dc === "dot-red") return "pill-red";
      if (dc === "dot-amber") return "pill-amber";
      return "";
    }

    function badgeToneClass(dotClassName) {
      if (dotClassName === "dot-green") return "pill-green";
      if (dotClassName === "dot-red") return "pill-red";
      if (dotClassName === "dot-amber") return "pill-amber";
      return "pill-gray";
    }

    function statusClass(status) {
      if (["OK", "Match", "MatchSTOP"].includes(status)) {
        return "status-ok";
      }
      if (["Translation mismatch", "Mismatch", "Error"].includes(status)) {
        return "status-bad";
      }
      return "status-na";
    }

    function toneClass(status) {
      if (["OK", "Match", "MatchSTOP"].includes(status)) {
        return "pill-green";
      }
      if (["Translation mismatch", "Mismatch", "Error"].includes(status)) {
        return "pill-red";
      }
      return "pill-blue";
    }

    function metadataStatusLabel(available, compatible = true) {
      if (!available) {
        return "missing";
      }
      if (!compatible) {
        return "incompatible";
      }
      return "available";
    }

    function wrapSequence(sequence, lineLength = 70) {
      const lines = [];
      for (let index = 0; index < sequence.length; index += lineLength) {
        lines.push(sequence.slice(index, index + lineLength));
      }
      return lines.join("\\n");
    }

    function buildExportHref(content) {
      return `data:text/plain;charset=utf-8,${encodeURIComponent(content)}`;
    }

    function buildRangeSelectionRows(analyses) {
      return analyses.map((analysis) => {
        const activeRange = getActiveRange(analysis);
        const construct = buildConstruct(analysis.entry, activeRange);
        const custom = primaryCustomRange(analysis);
        const automatic = automaticRange(analysis);
        const customDiff = customAutoDifference(analysis);
        const source = state.manualRanges[analysis.entry.id]
          ? "manual"
          : defaultRangeSourceLabel(analysis);
        return {
          id: analysis.entry.id,
          source,
          aaStart: activeRange.start,
          aaEnd: activeRange.end,
          aaRange: formatRange(activeRange),
          aaLength: activeRange.end - activeRange.start + 1,
          customStart: custom?.start ?? "",
          customEnd: custom?.end ?? "",
          cdsStart: construct.cdsRange?.start ?? "",
          cdsEnd: construct.cdsRange?.end ?? "",
          cdsRange: construct.cdsRange ? formatRange(construct.cdsRange) : "",
          cdsLength: construct.cds.length || "",
          constructStatus: construct.status,
          cdsStatus: analysis.validation.lengthStatus,
          translationStatus: analysis.validation.translationStatus,
          automaticLabel: analysis.suggestion.recommendedKey ?? analysis.entry.reference?.picked_range_name ?? "auto",
          automaticStart: automatic?.start ?? "",
          automaticEnd: automatic?.end ?? "",
          automaticRange: automatic ? formatRange(automatic) : "",
          customRange: custom ? formatRange(custom) : "",
          customAuto: customDiff ? (customDiff.different ? "different" : "match") : "",
        };
      });
    }

    function buildRangesExportTsv(rows) {
      const headers = [
        "id",
        "source",
        "aa_start",
        "aa_end",
        "aa_range",
        "aa_length",
        "cds_start",
        "cds_end",
        "cds_range",
        "cds_length",
        "construct_status",
        "cds_status",
        "translation_status",
        "automatic_label",
        "automatic_range",
        "custom_range",
        "custom_vs_automatic",
      ];
      const lines = [headers.join("\\t")];
      rows.forEach((row) => {
        lines.push([
          row.id,
          row.source,
          row.aaStart,
          row.aaEnd,
          row.aaRange,
          row.aaLength,
          row.cdsStart,
          row.cdsEnd,
          row.cdsRange,
          row.cdsLength,
          row.constructStatus,
          row.cdsStatus,
          row.translationStatus,
          row.automaticLabel,
          row.automaticRange,
          row.customRange,
          row.customAuto,
        ].join("\\t"));
      });
      return lines.join("\\n");
    }

    const BROWSER_LABEL_W = 124;
    const BROWSER_SVG_W = 1040;
    const BROWSER_TRACK_W = BROWSER_SVG_W - BROWSER_LABEL_W;

    function renderTrackViewer(entry, candidates, activeRange, customRanges = []) {
      const LABEL_W = BROWSER_LABEL_W;
      const SVG_W = BROWSER_SVG_W;
      const TRACK_W = BROWSER_TRACK_W;
      const GAP = 2;
      const mergedDomains = effectiveMergedDomains(entry);
      const hasIndividualDomains = entry.individualDomains.length > 0;
      const individualLaneHeight = 9;
      const individualLaneGap = 2;

      function layoutIndividualDomains(domains) {
        const laneEnds = [];
        return domains
          .slice()
          .sort((a, b) => a.start - b.start || a.end - b.end || a.label.localeCompare(b.label))
          .map((domain) => {
            let laneIndex = laneEnds.findIndex((laneEnd) => domain.start > laneEnd);
            if (laneIndex === -1) {
              laneIndex = laneEnds.length;
              laneEnds.push(domain.end);
            } else {
              laneEnds[laneIndex] = domain.end;
            }
            return { ...domain, laneIndex };
          });
      }

      const individualDomains = hasIndividualDomains ? layoutIndividualDomains(entry.individualDomains) : [];
      const individualLaneCount = individualDomains.length
        ? Math.max(...individualDomains.map((domain) => domain.laneIndex + 1))
        : 0;
      const domainRowHeight = hasIndividualDomains
        ? (
            6 + // top padding
            7 + // individual row label
            3 + // gap after individual label
            individualLaneCount * individualLaneHeight +
            Math.max(0, individualLaneCount - 1) * individualLaneGap +
            5 // bottom padding
          )
        : 24;
      const RANGES_LANE_H = 11;
      const RANGES_GAP = 3;
      const RANGES_PAD = 6;
      const totalRangeLaneCount = Math.max(1, candidates.length + customRanges.length);
      const rangesRowHeight =
        RANGES_PAD +
        totalRangeLaneCount * RANGES_LANE_H +
        Math.max(0, totalRangeLaneCount - 1) * RANGES_GAP +
        RANGES_PAD;

      const rows = [
        { label: "", height: 27, type: "ruler" },
        { label: "Domains", height: domainRowHeight, type: "domains" },
        { label: "Ranges", height: rangesRowHeight, type: "ranges" }
      ];

      function xPos(pos) {
        if (entry.proteinSequence.length <= 1) {
          return LABEL_W;
        }
        return LABEL_W + ((pos - 1) / (entry.proteinSequence.length - 1)) * TRACK_W;
      }

      const majorInterval =
        entry.proteinSequence.length > 400 ? 100 : entry.proteinSequence.length > 150 ? 50 : 25;
      const pxPerAA = TRACK_W / Math.max(1, entry.proteinSequence.length - 1);
      const showMinorTicks = pxPerAA >= 1.5; // per-AA ticks only when there's visible space

      const trackTops = [];
      let totalHeight = 0;
      rows.forEach((row) => {
        trackTops.push(totalHeight);
        totalHeight += row.height + GAP;
      });

      const selectedStartX = xPos(activeRange.start);
      const selectedEndX = xPos(activeRange.end);
      const svg = [];
      let domainClipCounter = 0;

      function pushLabeledDomainBox(kind, domain, x, y, width, height, color, opacity, fontSize, titleText = null) {
        const clipId = `${kind}-domain-clip-${domainClipCounter++}`;
        svg.push(
          `<g><title>${escapeHtml(titleText ?? `${domain.label} ${domain.start}-${domain.end}`)}</title>` +
            `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="${color}" opacity="${opacity}" rx="1.5" stroke="#000" stroke-width="1" />` +
            `<defs><clipPath id="${clipId}"><rect x="${x + 1}" y="${y + 1}" width="${Math.max(0, width - 2)}" height="${Math.max(0, height - 2)}" rx="1.5" /></clipPath></defs>` +
            `<text x="${x + 3}" y="${y + height / 2 + fontSize / 3 - 0.5}" font-size="${fontSize}" fill="#111" clip-path="url(#${clipId})">${escapeHtml(domain.label)}</text>` +
          `</g>`
        );
      }

      function pushCustomRangeBox(customRange, x, y, width, height) {
        const clipId = `custom-range-clip-${domainClipCounter++}`;
        const rangeLabel = `${customRange.label}: ${customRange.start}–${customRange.end}`;
        const titleText = `${customRange.label} (aa ${customRange.start}-${customRange.end})`;
        svg.push(
          `<g><title>${escapeHtml(titleText)}</title>` +
            `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="#dddddd" opacity="0.85" rx="1.5" stroke="#7a7a7a" stroke-width="1.5" stroke-dasharray="5,3" />` +
            `<defs><clipPath id="${clipId}"><rect x="${x + 1}" y="${y + 1}" width="${Math.max(0, width - 2)}" height="${Math.max(0, height - 2)}" rx="1.5" /></clipPath></defs>` +
            `<text x="${x + 3}" y="${y + height / 2 + 8 / 3 - 0.5}" font-size="8" fill="#333" clip-path="url(#${clipId})">${escapeHtml(rangeLabel)}</text>` +
          `</g>`
        );
      }

      function candidateHoverText(candidate) {
        return `${candidate.key}: ${candidate.description} (aa ${candidate.start}-${candidate.end})`;
      }

      const domainRowTop = trackTops[1];
      const rangesRowTop = trackTops[2];
      const individualLabelY = domainRowTop + 12;
      const individualTrackTop = domainRowTop + 16;
      const individualBlockHeight =
        individualLaneCount * individualLaneHeight +
        Math.max(0, individualLaneCount - 1) * individualLaneGap;
      const domainNoteY = domainRowTop + 16;

      rows.forEach((row, index) => {
        const y = trackTops[index];
        const bg = index % 2 === 0 ? "#ffffff" : "#f8f8f8";
        svg.push(
          `<rect x="0" y="${y}" width="${SVG_W}" height="${row.height}" fill="${bg}" />`,
          `<rect x="0" y="${y}" width="${LABEL_W}" height="${row.height}" fill="#f0f0f0" />`,
          row.type === "domains" && hasIndividualDomains
            ? `<text x="${LABEL_W - 4}" y="${y + 12}" text-anchor="end" font-size="9" fill="#555">${escapeHtml(row.label)}</text>`
            : row.type !== "ruler"
              ? `<text x="${LABEL_W - 4}" y="${y + row.height / 2 + 4}" text-anchor="end" font-size="9" fill="#555">` +
                  `${escapeHtml(row.label)}</text>`
              : "",
          row.type === "domains" && hasIndividualDomains
            ? `<text x="${LABEL_W - 4}" y="${individualLabelY}" text-anchor="end" font-size="7" fill="#666">individual</text>`
            : "",
          index > 0 ? `<line x1="0" y1="${y}" x2="${SVG_W}" y2="${y}" stroke="#d0d0d0" stroke-width="1" />` : "",
          `<line x1="${LABEL_W}" y1="${y}" x2="${LABEL_W}" y2="${y + row.height}" stroke="#aaa" stroke-width="1" />`
        );
      });

      // range bars: one lane per candidate (r1 / r2 / r3)
      candidates.forEach((candidate, idx) => {
        const color = candidateColors[candidate.key] ?? "#888";
        const laneY = rangesRowTop + RANGES_PAD + idx * (RANGES_LANE_H + RANGES_GAP);
        const x = xPos(candidate.start);
        const width = Math.max(4, xPos(candidate.end) - x);
        const rangeLabel = `${candidate.key}: ${candidate.start}–${candidate.end}`;
        const domainObj = { label: rangeLabel, start: candidate.start, end: candidate.end };
        pushLabeledDomainBox(
          "candidate",
          domainObj,
          x,
          laneY,
          width,
          RANGES_LANE_H,
          color,
          0.75,
          8,
          candidateHoverText(candidate),
        );
      });

      customRanges.forEach((customRange, idx) => {
        const laneY =
          rangesRowTop +
          RANGES_PAD +
          (candidates.length + idx) * (RANGES_LANE_H + RANGES_GAP);
        const x = xPos(customRange.start);
        const width = Math.max(4, xPos(customRange.end) - x);
        pushCustomRangeBox(customRange, x, laneY, width, RANGES_LANE_H);
      });

      svg.push(
        `<rect x="${selectedStartX}" y="0" width="${Math.max(0, selectedEndX - selectedStartX)}"` +
          ` height="${totalHeight}" fill="rgba(44,62,80,0.10)" />`
      );

      const rulerY = trackTops[0];
      const rulerBaseline = rulerY + rows[0].height - 2;
      svg.push(`<line x1="${LABEL_W}" y1="${rulerBaseline}" x2="${SVG_W}" y2="${rulerBaseline}" stroke="#bbb" stroke-width="1" />`);

      // minor ticks: every 1 AA (only when there's visible space)
      if (showMinorTicks) {
        for (let pos = 1; pos <= entry.proteinSequence.length; pos++) {
          if (pos % 10 === 0 || pos === 1 || pos === entry.proteinSequence.length) continue;
          const x = xPos(pos);
          svg.push(`<line x1="${x}" y1="${rulerBaseline - 2}" x2="${x}" y2="${rulerBaseline}" stroke="#ccc" stroke-width="0.5" />`);
        }
      }

      // medium ticks: every 10 AA (with label)
      for (let pos = 10; pos <= entry.proteinSequence.length; pos += 10) {
        if (pos % majorInterval === 0) continue; // major tick will cover this
        const x = xPos(pos);
        svg.push(
          `<line x1="${x}" y1="${rulerBaseline - 4}" x2="${x}" y2="${rulerBaseline}" stroke="#aaa" stroke-width="0.75" />`,
          `<text x="${x}" y="${rulerBaseline - 6}" text-anchor="middle" font-size="7" fill="#999">${pos}</text>`
        );
      }

      // major ticks: labeled at majorInterval (and first/last)
      const majorTicks = [1];
      for (let pos = majorInterval; pos < entry.proteinSequence.length; pos += majorInterval) majorTicks.push(pos);
      if (majorTicks[majorTicks.length - 1] !== entry.proteinSequence.length) majorTicks.push(entry.proteinSequence.length);
      majorTicks.forEach((pos) => {
        const x = xPos(pos);
        svg.push(
          `<line x1="${x}" y1="${rulerBaseline - 7}" x2="${x}" y2="${rulerBaseline}" stroke="#999" stroke-width="1" />`,
          `<text x="${x}" y="${rulerBaseline - 9}" text-anchor="middle" font-size="8" fill="#666">${pos}</text>`
        );
      });

      if (!mergedDomains.length && !entry.individualDomains.length) {
        svg.push(
          `<text x="${LABEL_W + 10}" y="${domainNoteY}" font-size="10" fill="#888">no domain BED ranges</text>`
        );
      }

      if (hasIndividualDomains) {
        individualDomains.forEach((domain, index) => {
          const x = xPos(domain.start);
          const width = Math.max(3, xPos(domain.end) - x);
          const color = trackPalette[index % trackPalette.length];
          const y = individualTrackTop + domain.laneIndex * (individualLaneHeight + individualLaneGap);
          pushLabeledDomainBox("individual", domain, x, y, width, individualLaneHeight, color, 0.9, 7);
        });
      }

      candidates.forEach((candidate) => {
        const color = candidateColors[candidate.key] ?? "#888";
        svg.push(
          `<line x1="${xPos(candidate.start)}" y1="0" x2="${xPos(candidate.start)}" y2="${totalHeight}"` +
            ` stroke="${color}" stroke-width="1.5" opacity="0.9" />`,
          `<line x1="${xPos(candidate.end)}" y1="0" x2="${xPos(candidate.end)}" y2="${totalHeight}"` +
            ` stroke="${color}" stroke-width="1.5" opacity="0.9" />`
        );
      });

      customRanges.forEach((customRange) => {
        const titleText = `${customRange.label} (aa ${customRange.start}-${customRange.end})`;
        const startX = xPos(customRange.start);
        const endX = xPos(customRange.end);
        svg.push(
          `<g><title>${escapeHtml(titleText)}</title><line x1="${startX}" y1="0" x2="${startX}" y2="${totalHeight}" stroke="#7a7a7a" stroke-width="1.5" stroke-dasharray="5,4" opacity="0.9" /></g>`,
          `<g><title>${escapeHtml(titleText)}</title><line x1="${endX}" y1="0" x2="${endX}" y2="${totalHeight}" stroke="#7a7a7a" stroke-width="1.5" stroke-dasharray="5,4" opacity="0.9" /></g>`
        );
      });

      svg.push(
        `<line x1="${selectedStartX}" y1="0" x2="${selectedStartX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="3" />`,
        `<line x1="${selectedEndX}" y1="0" x2="${selectedEndX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="3" />`,
        `<circle cx="${selectedStartX}" cy="8" r="3.5" fill="${selectedRangeColor}" />`,
        `<circle cx="${selectedEndX}" cy="8" r="3.5" fill="${selectedRangeColor}" />`
      );

      svg.push(
        `<line x1="${selectedStartX}" y1="0" x2="${selectedStartX}" y2="${totalHeight}"` +
          ` stroke="transparent" stroke-width="14" data-drag-handle="start" />`,
        `<line x1="${selectedEndX}" y1="0" x2="${selectedEndX}" y2="${totalHeight}"` +
          ` stroke="transparent" stroke-width="14" data-drag-handle="end" />`
      );

      return `
        <svg
          id="range-track-svg"
          class="track-svg"
          viewBox="0 0 ${SVG_W} ${totalHeight}"
          preserveAspectRatio="none"
          data-label-width="${LABEL_W}"
          data-svg-width="${SVG_W}"
          data-track-width="${TRACK_W}"
          data-protein-length="${entry.proteinSequence.length}"
        >
          ${svg.join("")}
        </svg>
      `;
    }

    function renderEvidenceBrowser(entry, activeRange) {
      const availableTracks = [];
      if (entry.evidence.structureDssp?.compatible) {
        availableTracks.push({
          type: "structure",
          label: "DSSP SS",
          coverage: entry.evidence.structureDssp.coverage,
          values: entry.evidence.structureDssp.values,
          height: 20
        });
      }
      if (entry.evidence.structureUniprot?.compatible) {
        availableTracks.push({
          type: "structure",
          label: "UniProt SS",
          coverage: entry.evidence.structureUniprot.coverage,
          values: entry.evidence.structureUniprot.values,
          height: 20
        });
      }
      if (entry.evidence.conservationFull?.compatible) {
        availableTracks.push({
          type: "numeric",
          label: "Cons full",
          tooltip: "Full-length orthogroup conservation from a fresh full-protein realignment; higher values mean stronger conservation.",
          coverage: entry.evidence.conservationFull.coverage,
          values: entry.evidence.conservationFull.values,
          height: 34,
          fill: "rgba(58,133,200,0.28)",
          stroke: "#2868a0"
        });
      }
      if (entry.evidence.conservation?.compatible) {
        availableTracks.push({
          type: "numeric",
          label: "Cons DBD",
          coverage: entry.evidence.conservation.coverage,
          values: entry.evidence.conservation.values,
          height: 34,
          fill: "rgba(39,174,96,0.28)",
          stroke: "#1f7a45"
        });
      }
      if (entry.evidence.iupred?.compatible) {
        availableTracks.push({
          type: "numeric",
          label: "IUPred",
          coverage: entry.evidence.iupred.coverage,
          values: entry.evidence.iupred.values,
          height: 34,
          fill: "rgba(110,110,110,0.20)",
          stroke: "#666666",
          maxValue: entry.evidence.iupred.maxValue ?? 1
        });
      }
      if (entry.evidence.plddt?.compatible) {
        availableTracks.push({
          type: "numeric",
          label: "pLDDT",
          coverage: entry.evidence.plddt.coverage,
          values: entry.evidence.plddt.values,
          height: 34,
          fill: "rgba(215,150,0,0.25)",
          stroke: "#c47d00",
          maxValue: entry.evidence.plddt.maxValue ?? 100
        });
      }

      if (!availableTracks.length) {
        return `<div class="browser-empty">No compatible evidence tracks are available for this protein.</div>`;
      }

      const LABEL_W = BROWSER_LABEL_W;
      const SVG_W = BROWSER_SVG_W;
      const TRACK_W = BROWSER_TRACK_W;
      const GAP = 2;
      const residueWidth =
        entry.proteinSequence.length > 1
          ? TRACK_W / (entry.proteinSequence.length - 1)
          : TRACK_W;

      function xPos(pos) {
        if (entry.proteinSequence.length <= 1) {
          return LABEL_W;
        }
        return LABEL_W + ((pos - 1) / (entry.proteinSequence.length - 1)) * TRACK_W;
      }

      const rows = [{ type: "ruler", label: "", height: 22 }, ...availableTracks];

      const trackTops = [];
      let totalHeight = 0;
      rows.forEach((row) => {
        trackTops.push(totalHeight);
        totalHeight += row.height + GAP;
      });

      const selectedStartX = xPos(activeRange.start);
      const selectedEndX = xPos(activeRange.end);
      const tickInterval =
        entry.proteinSequence.length > 400 ? 100 : entry.proteinSequence.length > 150 ? 50 : 25;
      const ticks = [1];
      for (let pos = tickInterval; pos < entry.proteinSequence.length; pos += tickInterval) {
        ticks.push(pos);
      }
      if (ticks[ticks.length - 1] !== entry.proteinSequence.length) {
        ticks.push(entry.proteinSequence.length);
      }

      const svg = [];
      rows.forEach((row, index) => {
        const y = trackTops[index];
        const bg = index % 2 === 0 ? "#ffffff" : "#f8f8f8";
        svg.push(
          `<rect x="0" y="${y}" width="${SVG_W}" height="${row.height}" fill="${bg}" />`,
          `<rect x="0" y="${y}" width="${LABEL_W}" height="${row.height}" fill="#efefef" />`,
          row.type !== "ruler"
            ? (
                row.tooltip
                  ? `<g><title>${escapeHtml(row.tooltip)}</title><text x="${LABEL_W - 6}" y="${y + row.height / 2 - 1}" text-anchor="end" font-size="9" fill="#444" font-weight="700">${escapeHtml(row.label)}</text></g>`
                  : `<text x="${LABEL_W - 6}" y="${y + row.height / 2 - 1}" text-anchor="end" font-size="9" fill="#444" font-weight="700">${escapeHtml(row.label)}</text>`
              )
            : "",
          row.type !== "ruler"
            ? `<text x="${LABEL_W - 6}" y="${y + row.height / 2 + 6}" text-anchor="end" font-size="8" fill="#777">${Math.round(row.coverage * 100)}%</text>`
            : "",
          index > 0 ? `<line x1="0" y1="${y}" x2="${SVG_W}" y2="${y}" stroke="#d5d5d5" stroke-width="1" />` : "",
          `<line x1="${LABEL_W}" y1="${y}" x2="${LABEL_W}" y2="${y + row.height}" stroke="#aaa" stroke-width="1" />`
        );
      });

      svg.push(
        `<rect x="${selectedStartX}" y="0" width="${Math.max(0, selectedEndX - selectedStartX)}"` +
          ` height="${totalHeight}" fill="rgba(44,62,80,0.14)" />`
      );

      const rulerY = trackTops[0];
      const rulerBaseline = rulerY + rows[0].height - 3;
      svg.push(
        `<line x1="${LABEL_W}" y1="${rulerBaseline}" x2="${SVG_W}" y2="${rulerBaseline}" stroke="#bcbcbc" stroke-width="1" />`
      );
      ticks.forEach((pos) => {
        const x = xPos(pos);
        svg.push(
          `<line x1="${x}" y1="${rulerBaseline - 4}" x2="${x}" y2="${rulerBaseline}" stroke="#999" stroke-width="1" />`,
          `<text x="${x}" y="${rulerBaseline - 6}" text-anchor="middle" font-size="8" fill="#666">${pos}</text>`
        );
      });

      availableTracks.forEach((track, index) => {
        const y = trackTops[index + 1];
        if (track.type === "structure") {
          buildRuns(track.values).forEach((run) => {
            if (run.value === null) {
              return;
            }
            const x = xPos(run.start + 1);
            const width = Math.max(2, (run.end - run.start + 1) * residueWidth);
            const color = ssColors[run.value] ?? "#bbbbbb";
            const tooltip = `${track.label}: ${structureTypeName(run.value)} (${run.value}), aa ${run.start + 1}-${run.end + 1}`;
            svg.push(
              `<rect x="${x}" y="${y + 3}" width="${width}" height="${track.height - 6}" fill="${color}" rx="1">` +
                `<title>${escapeHtml(tooltip)}</title>` +
              `</rect>`
            );
          });
        } else {
          const trackTop = y + 3;
          const trackHeight = track.height - 6;
          const maxValue = track.maxValue != null
            ? track.maxValue
            : Math.max(
                1,
                ...track.values.filter((value) => value !== null && value !== undefined).map(Number)
              );
          const areaPath = buildAreaPath(track.values, xPos, trackTop, trackHeight, maxValue);
          // y-axis grid lines at 25/50/75% and labels at top/mid/bottom
          const ySteps = [0, 0.25, 0.5, 0.75, 1.0];
          ySteps.forEach((frac) => {
            const lineY = trackTop + trackHeight * (1 - frac);
            svg.push(
              `<line x1="${LABEL_W}" y1="${lineY}" x2="${SVG_W}" y2="${lineY}" stroke="#e0e0e0" stroke-width="0.5" stroke-dasharray="3,3" />`
            );
          });
          // y-axis value labels: top, mid, bottom
          const yAxisLabels = [
            { frac: 1.0, val: maxValue },
            { frac: 0.5, val: maxValue / 2 },
            { frac: 0.0, val: 0 }
          ];
          yAxisLabels.forEach(({ frac, val }) => {
            const labelY = trackTop + trackHeight * (1 - frac);
            const text = Number.isInteger(val) ? String(val) : val.toFixed(1);
            svg.push(
              `<text x="${LABEL_W - 3}" y="${labelY + 3}" text-anchor="end" font-size="7" fill="#999">${text}</text>`
            );
          });
          if (areaPath) {
            svg.push(
              `<path d="${areaPath}" fill="${track.fill}" stroke="${track.stroke}" stroke-width="1.2" />`
            );
          }
          svg.push(
            `<line x1="${LABEL_W}" y1="${y + track.height - 4}" x2="${SVG_W}" y2="${y + track.height - 4}" stroke="#d0d0d0" stroke-width="1" />`
          );
        }
      });

      svg.push(
        `<line x1="${selectedStartX}" y1="0" x2="${selectedStartX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="2" />`,
        `<line x1="${selectedEndX}" y1="0" x2="${selectedEndX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="2" />`
      );

      return `
        <div class="browser-wrap">
          <svg class="browser-svg" viewBox="0 0 ${SVG_W} ${totalHeight}" preserveAspectRatio="none">
            ${svg.join("")}
          </svg>
          <div class="browser-legend">
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-helix"></span>helix / structured helix-like</span>
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-strand"></span>strand / sheet-like</span>
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-coil"></span>coil / turn classes</span>
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-numeric"></span>numeric tracks use their own scale (conservation, IUPred 0–1, pLDDT 0–100)</span>
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-range"></span>selected construct span</span>
          </div>
        </div>
      `;
    }

    function renderStructurePanel(entry, activeRange, docked = false) {
      const model = entry.evidence.structureModel;
      if (!model?.text) {
        return "";
      }
      const selection = getStructureSelection(entry, activeRange);

      return `
        <div class="structure-pane-wrap ${docked ? "structure-pane-wrap-docked" : ""}">
          ${
            docked
              ? ""
              : `<div class="structure-pane-meta">
                  <span class="metric-chip">AA ${escapeHtml(formatRange(activeRange))}</span>
                  ${
                    selection.hasMapping && selection.proteinRange
                      ? `<span class="metric-chip">model aa ${escapeHtml(formatRange(selection.proteinRange))}</span>`
                      : `<span class="metric-chip">model mapping unavailable</span>`
                  }
                  ${
                    selection.mappingSource
                      ? `<span class="metric-chip">mapped by ${escapeHtml(selection.mappingSource)}</span>`
                      : ""
                  }
                  ${
                    selection.hasOverlap && selection.overlapRange
                      ? `<span class="metric-chip">shown ${escapeHtml(formatRange(selection.overlapRange))}</span>`
                      : selection.hasMapping
                        ? `<span class="metric-chip">no overlap with selected range</span>`
                        : ""
                  }
                  <span class="metric-chip">${escapeHtml(model.source)}</span>
                </div>`
          }
          <div class="structure-panel ${docked ? "structure-panel-docked" : ""}">
            <div
              id="structure-viewer-${sanitizeDomId(entry.id)}"
              class="structure-viewer"
              data-structure-viewer
            ></div>
          </div>
        </div>
      `;
    }

    function _applyStructureColors(viewer, entry, activeRange, withZoom) {
      viewer.setStyle({}, { cartoon: { color: "#d0d0d0" } });
      const selection = getStructureSelection(entry, activeRange);
      function rainbowColor(fraction) {
        const stops = [
          [217, 55, 55],
          [235, 140, 30],
          [232, 196, 54],
          [65, 170, 95],
          [52, 120, 210],
          [130, 80, 190],
        ];
        const clampedFraction = clamp(fraction, 0, 1);
        const scaled = clampedFraction * (stops.length - 1);
        const index = Math.min(stops.length - 2, Math.floor(scaled));
        const localFraction = scaled - index;
        const start = stops[index];
        const end = stops[index + 1];
        const rgb = start.map((channel, channelIndex) =>
          Math.round(channel + (end[channelIndex] - channel) * localFraction)
        );
        return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
      }

      if (selection.hasOverlap && selection.residueRange) {
        const residueNumbers = Array.isArray(entry.evidence.structureModel?.mapping?.residueNumbers)
          ? entry.evidence.structureModel.mapping.residueNumbers
          : [];
        const residueIndexByNumber = new Map(
          residueNumbers.map((residueNumber, index) => [Number(residueNumber), index])
        );
        const proteinStart = Number(entry.evidence.structureModel?.mapping?.proteinStart) || 1;
        const proteinEnd = Number(entry.evidence.structureModel?.mapping?.proteinEnd) || proteinStart;
        const denominator = Math.max(1, proteinEnd - proteinStart);
        viewer.setStyle(
          { resi: `${selection.residueRange.start}-${selection.residueRange.end}` },
          {
            cartoon: {
              colorfunc: (atom) => {
                const residueIndex = residueIndexByNumber.get(Number(atom.resi));
                if (residueIndex === undefined) {
                  return rainbowColor(0);
                }
                return rainbowColor(residueIndex / denominator);
              }
            }
          }
        );
        if (withZoom) {
          viewer.zoomTo({ resi: `${selection.residueRange.start}-${selection.residueRange.end}` });
        }
      } else if (withZoom) {
        viewer.zoomTo();
      }
    }

    function stabilizeStructureViewer(viewer) {
      const refresh = () => {
        viewer.resize();
        viewer.render();
      };
      refresh();
      window.requestAnimationFrame(() => {
        refresh();
        window.setTimeout(refresh, 40);
      });
    }

    function initializeStructureViewer(entry, activeRange, docked = false) {
      const layoutMode = docked ? "docked" : "stacked";
      const initToken = ++_structureInitToken;

      function tryMount(attempt = 0) {
        if (initToken !== _structureInitToken) {
          return;
        }

        const newSlot = detailPanelEl.querySelector("[data-structure-viewer]");
        const model = entry.evidence.structureModel;
        if (!newSlot || !model?.text) {
          _structureCache = { entryId: null, viewer: null, containerEl: null, layoutMode: null };
          return;
        }

        const rect = newSlot.getBoundingClientRect();
        const hasRealViewport = rect.width >= 80 && rect.height >= 120;
        if (!hasRealViewport && attempt < 8) {
          window.setTimeout(() => tryMount(attempt + 1), 30);
          return;
        }

        if (
          _structureCache.entryId === entry.id &&
          _structureCache.layoutMode === layoutMode &&
          _structureCache.viewer &&
          _structureCache.containerEl
        ) {
          newSlot.replaceWith(_structureCache.containerEl);
          _applyStructureColors(_structureCache.viewer, entry, activeRange, false);
          stabilizeStructureViewer(_structureCache.viewer);
          return;
        }

        _structureCache = { entryId: null, viewer: null, containerEl: null, layoutMode: null };

        if (!window.$3Dmol || typeof window.$3Dmol.createViewer !== "function") {
          newSlot.innerHTML = `
            <div class="structure-fallback">
              The 3D viewer library did not load. Open the report with internet access to fetch 3Dmol.js, or vendor the script locally for offline viewing.
            </div>
          `;
          return;
        }

        try {
          const viewer = window.$3Dmol.createViewer(newSlot, { backgroundColor: "white" });
          viewer.addModel(model.text, model.format || "pdb");
          _applyStructureColors(viewer, entry, activeRange, true);
          stabilizeStructureViewer(viewer);
          window.setTimeout(() => {
            if (initToken !== _structureInitToken) {
              return;
            }
            stabilizeStructureViewer(viewer);
          }, 30);
          _structureCache = { entryId: entry.id, viewer, containerEl: newSlot, layoutMode };
        } catch (error) {
          const message = error instanceof Error ? error.message : String(error);
          newSlot.innerHTML = `
            <div class="structure-fallback">
              Failed to render the local PDB model: ${escapeHtml(message)}
            </div>
          `;
        }
      }

      window.setTimeout(() => tryMount(), 0);
    }

    function hasCdsIssue(analysis, construct) {
      return (
        analysis.validation.translationStatus !== "Match" ||
        !["Match", "MatchSTOP"].includes(analysis.validation.lengthStatus) ||
        construct.status !== "OK"
      );
    }

    function renderComparisonBlocks(expected, observed, cdsSequence, expectedLabel, observedLabel) {
      const chunkSize = 40;
      const maxLen = Math.max(expected.length, observed.length);
      const blocks = [];

      function renderAaChar(char, mismatch) {
        const display = char || "·";
        const cls = mismatch
          ? (char ? "issue-char issue-char-mismatch" : "issue-char issue-char-missing")
          : "issue-char";
        return `<span class="${cls}">${escapeHtml(display)}</span>`;
      }

      function renderDiffChar(mismatch) {
        return `<span class="issue-diff-char ${mismatch ? "issue-diff-mismatch" : "issue-diff-match"}">${mismatch ? "^" : "|"}</span>`;
      }

      for (let chunkStart = 0; chunkStart < maxLen; chunkStart += chunkSize) {
        const chunkEnd = Math.min(maxLen, chunkStart + chunkSize);
        let expectedHtml = "";
        let observedHtml = "";
        let diffHtml = "";
        let cdsHtml = "";

        for (let index = chunkStart; index < chunkEnd; index += 1) {
          const expectedChar = expected[index] ?? "";
          const observedChar = observed[index] ?? "";
          const mismatch = expectedChar !== observedChar;
          expectedHtml += renderAaChar(expectedChar, mismatch);
          observedHtml += renderAaChar(observedChar, mismatch);
          diffHtml += renderDiffChar(mismatch);

          if (cdsSequence) {
            const codon = cdsSequence.slice(index * 3, index * 3 + 3);
            if (codon) {
              cdsHtml += `<span class="${mismatch ? "issue-codon issue-codon-mismatch" : "issue-codon"}">${escapeHtml(codon)}</span>`;
            }
          }
        }

        blocks.push(`
          <div class="issue-compare-block">
            <div class="issue-compare-row">
              <span class="issue-row-label">${escapeHtml(expectedLabel)}</span>
              <span class="issue-row-start">${chunkStart + 1}</span>
              <span class="issue-row-seq">${expectedHtml}</span>
            </div>
            <div class="issue-compare-row">
              <span class="issue-row-label">${escapeHtml(observedLabel)}</span>
              <span class="issue-row-start">${chunkStart + 1}</span>
              <span class="issue-row-seq">${observedHtml}</span>
            </div>
            <div class="issue-compare-row">
              <span class="issue-row-label">diff</span>
              <span class="issue-row-start"></span>
              <span class="issue-row-seq issue-row-seq-diff">${diffHtml}</span>
            </div>
            ${
              cdsHtml
                ? `<div class="issue-compare-row">
                    <span class="issue-row-label">CDS</span>
                    <span class="issue-row-start">${chunkStart * 3 + 1}</span>
                    <span class="issue-row-cds">${cdsHtml}</span>
                  </div>`
                : ""
            }
          </div>
        `);
      }

      return `<div class="issue-compare">${blocks.join("")}</div>`;
    }

    function renderCdsIssuePanel(entry, analysis, activeRange, construct) {
      const issueSections = [];
      const requestedCdsRange = toCdsRange(activeRange);
      const fullTranslation = translateDna(entry.cdsSequence).replace(/\*$/, "");
      const lengthProblem = !["Match", "MatchSTOP"].includes(analysis.validation.lengthStatus);
      const constructOverflow = construct.status === "Range exceeds CDS length";

      if (lengthProblem) {
        const expectedNt = entry.proteinSequence.length * 3;
        const actualNt = entry.cdsSequence.length;
        issueSections.push(`
          <div class="issue-alert issue-alert-warning">
            Full CDS length does not match the protein record: protein implies ${expectedNt} nt, but the CDS file contains ${actualNt} nt.
          </div>
        `);
      }

      if (analysis.validation.translationStatus !== "Match") {
        const detailParts = [];
        if (analysis.validation.firstMismatchAa !== null) {
          detailParts.push(`first mismatch at aa ${analysis.validation.firstMismatchAa}`);
        }
        if (analysis.validation.lengthDelta !== 0) {
          detailParts.push(`CDS translates to ${analysis.validation.cdsPeptideLength} aa (${analysis.validation.lengthDelta > 0 ? "+" : ""}${analysis.validation.lengthDelta} vs protein)`);
        }
        issueSections.push(`
          <div class="issue-subsection">
            <div class="controls-section-label">Full protein vs CDS translation</div>
            <div class="issue-subsection-note">
              The full CDS does not translate exactly to the protein record${detailParts.length ? `; ${escapeHtml(detailParts.join(", "))}` : ""}.
            </div>
            ${renderComparisonBlocks(entry.proteinSequence, fullTranslation, entry.cdsSequence, "protein", "translated")}
          </div>
        `);
      }

      if (construct.status === "Translation mismatch") {
        issueSections.push(`
          <div class="issue-subsection">
            <div class="controls-section-label">Selected construct vs selected CDS translation</div>
            <div class="issue-subsection-note">
              The currently selected range ${escapeHtml(formatRange(activeRange))} does not translate cleanly from the selected CDS span ${escapeHtml(formatRange(construct.cdsRange))}.
            </div>
            ${renderComparisonBlocks(construct.peptide, construct.translation, construct.cds, "selected AA", "translated")}
          </div>
        `);
      }

      if (constructOverflow) {
        const cdsStart = (activeRange.start - 1) * 3;
        const availableTail = entry.cdsSequence.slice(Math.min(cdsStart, entry.cdsSequence.length));
        const missingNt = Math.max(0, requestedCdsRange.end - entry.cdsSequence.length);
        issueSections.push(`
          <div class="issue-subsection">
            <div class="controls-section-label">Selected range exceeds CDS length</div>
            <div class="issue-alert issue-alert-danger">
              The selected AA range ${escapeHtml(formatRange(activeRange))} requires CDS coordinates ${escapeHtml(formatRange(requestedCdsRange))}, but the CDS ends at nt ${entry.cdsSequence.length}.
            </div>
            <div class="issue-subsection-note">
              The highlighted CDS tail below is all that exists from the requested start; ${missingNt} nt are missing beyond the CDS end.
            </div>
            <div class="seq-panel">
              <div class="seq-track-label">Selected peptide</div>
              <pre class="seq-pre">${escapeHtml(construct.peptide || entry.proteinSequence.slice(activeRange.start - 1, activeRange.end))}</pre>
              <div class="seq-track-label">Available CDS from requested start</div>
              <pre class="seq-pre"><mark class="seq-mark">${escapeHtml(availableTail)}</mark><span class="issue-inline-missing">[missing ${missingNt} nt beyond CDS end]</span></pre>
            </div>
          </div>
        `);
      }

      if (!issueSections.length) {
        return "";
      }

      return `
        <div class="issue-panel-body">
          <div class="issue-legend">
            <span class="issue-legend-item"><span class="issue-legend-swatch issue-legend-swatch-mismatch"></span>mismatch</span>
            <span class="issue-legend-item"><span class="issue-legend-swatch issue-legend-swatch-missing"></span>missing residue / truncated end</span>
            <span class="issue-legend-item">diff row: <strong>^</strong> mismatch, <strong>|</strong> match</span>
            <span class="issue-legend-item">missing residues are shown as <strong>·</strong></span>
          </div>
          ${issueSections.join("")}
        </div>
      `;
    }

    function renderToolbar(selectedAnalysis, analyses) {
      if (!selectedAnalysis) {
        toolbarEl.innerHTML = `
          <div class="app-toolbar-title">
            <strong>Construct Design HTML Report</strong>
            <span>No proteins matched the current filters.</span>
          </div>
        `;
        return;
      }

      const activeRange = getActiveRange(selectedAnalysis);
      const construct = buildConstruct(selectedAnalysis.entry, activeRange);
      const selectedExportHref =
        construct.status === "OK"
          ? buildExportHref(
              `>${selectedAnalysis.entry.id}|aa=${formatRange(activeRange)}\\n${construct.cds}`
            )
          : null;

      const recommendedBatchFasta = analyses
        .flatMap((analysis) => {
          const activeRange = getActiveRange(analysis);
          if (!activeRange) {
            return [];
          }

          const itemConstruct = buildConstruct(analysis.entry, activeRange);
          if (itemConstruct.status !== "OK") {
            return [];
          }

          const rangeLabel = state.manualRanges[analysis.entry.id]
            ? "manual"
            : defaultRangeSourceLabel(analysis);
          return [
            `>${analysis.entry.id}|${rangeLabel}|aa=${formatRange(activeRange)}\\n${itemConstruct.cds}`
          ];
        })
        .join("\\n");

      toolbarEl.innerHTML = `
        <div class="app-toolbar-title">
          <strong>Construct Design HTML Report</strong>
          <span>Required inputs: protein FASTA, CDS FASTA, optional domain BED</span>
          <span>Records: ${escapeHtml(datasetSummary.display)}</span>
          <span>Inputs: ${escapeHtml(payload.inputSummary)} · generated ${escapeHtml(payload.generatedAt)}</span>
        </div>
        <div class="app-toolbar-actions">
          ${
            selectedExportHref
              ? `<a class="action-button" href="${selectedExportHref}" download="${escapeHtml(selectedAnalysis.entry.id)}.construct.fasta">Export selected FASTA</a>`
              : `<span class="action-button action-button-muted">Not export-ready</span>`
          }
          <a
            class="action-button action-button-secondary"
            href="${buildExportHref(recommendedBatchFasta)}"
            download="recommended_constructs.fasta"
          >
            Export batch
          </a>
        </div>
      `;
    }

    function renderBatch(visibleAnalyses, selectedId) {
      batchPanelEl.innerHTML = `
        <div class="batch-panel-header">
          <h2>Proteins (${visibleAnalyses.length})</h2>
          <p>Compact navigator for selecting a protein to inspect in detail.</p>
        </div>
        <div class="batch-toolbar">
          <input id="search-input" value="${escapeHtml(state.search)}" placeholder="Search sequence ID or domain…" aria-label="Search">
          <select id="filter-input" aria-label="Filter">
            <option value="all"${state.filter === "all" ? " selected" : ""}>All</option>
            <option value="good"${state.filter === "good" ? " selected" : ""}>Good</option>
            <option value="range_mismatch"${state.filter === "range_mismatch" ? " selected" : ""}>Range mismatch</option>
            <option value="cds_error"${state.filter === "cds_error" ? " selected" : ""}>CDS error</option>
          </select>
          <select id="sort-input" aria-label="Sort by">
            <option value="id"${state.sort === "id" ? " selected" : ""}>Sort: ID</option>
            <option value="evidence"${state.sort === "evidence" ? " selected" : ""}>Sort: evidence</option>
            <option value="structure"${state.sort === "structure" ? " selected" : ""}>Sort: 3D+SS</option>
            <option value="conservation"${state.sort === "conservation" ? " selected" : ""}>Sort: conservation</option>
            <option value="disorder"${state.sort === "disorder" ? " selected" : ""}>Sort: IUPred</option>
            <option value="plddt"${state.sort === "plddt" ? " selected" : ""}>Sort: pLDDT</option>
          </select>
        </div>
        <div class="review-list">
          ${
            visibleAnalyses.length
              ? visibleAnalyses
                  .map((analysis) => {
                    const active = analysis.entry.id === selectedId;
                    const activeRange = getActiveRange(analysis);
                    const evidenceCardsForRow = evidenceCards(analysis.entry, analysis.validation);
                    const customDiff = customAutoDifference(analysis);
                    const availableEvidenceCount = evidenceCardsForRow.filter((card) => card.cls !== "dot-gray").length;
                    const issueCount = evidenceCardsForRow.filter((card) => card.cls === "dot-red").length;
                    const domainSummary =
                      analysis.entry.individualDomains.map((d) => d.label).join(", ") ||
                      analysis.entry.mergedDomains.map((d) => d.label).join(", ") ||
                      "no domains";
                    const rangeSource = state.manualRanges[analysis.entry.id]
                      ? "manual"
                      : defaultRangeSourceLabel(analysis);
                    const rangeStatus = analysisReviewStatus(analysis);
                    const groupStatus = analysisGroup(analysis);
                    const structureStatus = analysis.entry.reference?.status_structure ?? "";
                    const cdsProblem = analysisHasCdsProteinProblem(analysis);

                    return `
                      <button type="button" class="review-row ${active ? "review-row-active" : ""}" data-select-id="${escapeHtml(analysis.entry.id)}" data-range-status="${escapeHtml(groupStatus === "Good" ? rangeStatus : groupStatus)}">
                        <div class="review-row-top">
                          <div class="review-row-title">
                            <strong>${escapeHtml(analysis.entry.id)}</strong>
                            <div class="review-row-meta">
                              <span class="review-mini-pill">${analysis.entry.proteinSequence.length} aa</span>
                              <span class="review-mini-pill ${state.manualRanges[analysis.entry.id] ? "pill-blue" : ""}">${escapeHtml(rangeSource)}: ${escapeHtml(formatRange(activeRange))}</span>
                              ${rangeStatus ? `<span class="review-mini-pill ${dotToTone(rangeStatus)}">range: ${escapeHtml(rangeStatus)}</span>` : ""}
                            </div>
                          </div>
                        </div>
                        <div class="review-row-id">${escapeHtml(domainSummary)}</div>
                        <div class="review-row-summary">
                          <span class="review-mini-pill ${
                            groupStatus === "Good" ? "pill-green" : groupStatus === "CDS error" ? "pill-red" : "pill-amber"
                          }">${escapeHtml(groupStatus)}</span>
                          <span class="review-mini-pill">evidence: ${availableEvidenceCount}/${evidenceCardsForRow.length}</span>
                          ${cdsProblem ? `<span class="review-mini-pill pill-red">CDS mismatch</span>` : ""}
                          ${
                            customDiff?.different
                              ? `<span class="review-mini-pill pill-amber" title="${escapeHtml(`custom ${formatRange(customDiff.custom)} differs from automatic ${customDiff.automaticLabel} ${formatRange(customDiff.automatic)}`)}">custom != ${escapeHtml(customDiff.automaticLabel)}</span>`
                              : ""
                          }
                          ${issueCount ? `<span class="review-mini-pill pill-red">${issueCount} issue${issueCount === 1 ? "" : "s"}</span>` : ""}
                          ${structureStatus ? `<span class="review-mini-pill ${dotToTone(structureStatus)}">structure: ${escapeHtml(structureStatus)}</span>` : ""}
                        </div>
                      </button>
                    `;
                  })
                  .join("")
              : `<div class="empty-state">No proteins matched the current filters.</div>`
          }
        </div>
      `;

      batchPanelEl.querySelector("#search-input")?.addEventListener("input", (event) => {
        state.search = event.target.value;
        debounceRender();
      });

      batchPanelEl.querySelector("#filter-input")?.addEventListener("change", (event) => {
        state.filter = event.target.value;
        render();
      });

      batchPanelEl.querySelector("#sort-input")?.addEventListener("change", (event) => {
        state.sort = event.target.value;
        render();
      });

      batchPanelEl.querySelectorAll("[data-select-id]").forEach((button) => {
        button.addEventListener("click", () => {
          state.selectedId = button.getAttribute("data-select-id");
          render();
        });
      });
    }

    function renderDetail(analysis) {
      if (!analysis) {
        detailPanelEl.innerHTML = `<div class="empty-state">No protein is currently selected.</div>`;
        return;
      }

      const entry = analysis.entry;
      const activeRange = getActiveRange(analysis);
      const construct = buildConstruct(entry, activeRange);
      const candidateRanges = Object.values(analysis.suggestion.candidates).filter(Boolean);
      const fastaPreview =
        construct.status === "OK"
          ? `>${entry.id}|aa=${formatRange(activeRange)}|cds=${formatRange(construct.cdsRange)}\\n${wrapSequence(construct.cds)}`
          : "";
      const aaLen = activeRange.end - activeRange.start + 1;
      const cdsLen = construct.cds.length;
      const detailEvidenceCards = evidenceCards(entry, analysis.validation);
      const hasStructureModel = Boolean(entry.evidence.structureModel?.text);
      const showCdsIssuePanel = hasCdsIssue(analysis, construct);
      const structureToggleLabel = state.structureDockRight ? "3D below" : "3D right";
      let sectionIndex = 1;
      const coordinatesSectionIndex = sectionIndex++;
      const evidenceSectionIndex = sectionIndex++;
      const structureSectionIndex = hasStructureModel ? sectionIndex++ : null;
      const cdsIssueSectionIndex = showCdsIssuePanel ? sectionIndex++ : null;
      const constructSummarySectionIndex = sectionIndex++;
      const coordinatesSectionMarkup = `
        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">${coordinatesSectionIndex}. Coordinates and Sequence</div>
              <div class="detail-section-copy">Protein coordinates with direct range editing controls.</div>
            </div>
            <div class="section-chip-row">
              <span class="metric-chip">selected: ${escapeHtml(formatRange(activeRange))}</span>
              <span class="metric-chip">${aaLen} aa</span>
              <span class="metric-chip">${cdsLen} nt</span>
              ${
                hasStructureModel
                  ? `<button type="button" class="coordinates-extra-toggle" data-toggle-structure-dock>${structureToggleLabel}</button>`
                  : ""
              }
              <span class="metric-chip">expand / collapse</span>
            </div>
          </summary>

          <div class="track-viewer-wrap">
            ${renderTrackViewer(entry, candidateRanges, activeRange, analysis.customRanges)}
            <div class="track-drag-note">Drag the black range boundaries directly on the coordinate plot.</div>
          </div>
          <div class="controls-row">
            <div>
              <div class="controls-section-label">Suggested ranges</div>
              <div class="candidate-buttons">
                ${
                  candidateRanges.length
                    ? candidateRanges
                        .map((candidate) => {
                          const isActive =
                            candidate.start === activeRange.start && candidate.end === activeRange.end;
                          return `
                            <button
                              type="button"
                              class="candidate-btn ${isActive ? "candidate-btn-active" : ""}"
                              data-candidate-range="${candidate.key}"
                              title="${escapeHtml(candidate.description)}"
                            >
                              <span class="candidate-btn-dot" style="background: ${candidateColors[candidate.key] ?? "#888"}"></span>
                              ${escapeHtml(candidate.label)} ${candidate.start}-${candidate.end}
                            </button>
                          `;
                        })
                        .join("")
                    : `<span class="metric-chip">No candidate ranges available</span>`
                }
              </div>
            </div>

            <div>
              <div class="controls-section-label">Manual range</div>
              <div class="manual-range-row">
                <label>Start
                  <input type="number" id="range-start-input" min="1" max="${entry.proteinSequence.length}" value="${activeRange.start}">
                </label>
                <label>End
                  <input type="number" id="range-end-input" min="${activeRange.start}" max="${entry.proteinSequence.length}" value="${activeRange.end}">
                </label>
              </div>
              <div class="range-status-line">
                construct status: <span class="${statusClass(construct.status)}">${escapeHtml(construct.status)}</span>
              </div>
            </div>
          </div>
        </details>
      `;
      const structureDockedMarkup = hasStructureModel
        ? `
          <div class="structure-docked-host">
            ${renderStructurePanel(entry, activeRange, true)}
          </div>
        `
        : "";
      const structureSectionMarkup = hasStructureModel
        ? `
          <details class="detail-section" open>
            <summary class="detail-section-header detail-section-summary">
              <div>
                <div class="detail-section-index">${structureSectionIndex}. Structure</div>
                <div class="detail-section-copy">Local protein structure with selected residues colored by their fixed N-to-C position in the full model, from N-terminal red through a C-terminal rainbow.</div>
              </div>
              <div class="metric-chip">expand / collapse</div>
            </summary>
            <div class="metadata-body">
              ${renderStructurePanel(entry, activeRange, state.structureDockRight)}
            </div>
          </details>
        `
        : "";
      const evidenceSectionMarkup = `
        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">${evidenceSectionIndex}. Evidence Tracks</div>
              <div class="detail-section-copy">Mapped browser tracks aligned to protein coordinates.</div>
            </div>
            <div class="metric-chip">expand / collapse</div>
          </summary>
          <div class="metadata-body">
            ${renderEvidenceBrowser(entry, activeRange)}
            <div class="metadata-note">
              Only tracks that were successfully loaded and mapped onto the current protein sequence are shown here.
            </div>
            ${
              entry.individualDomains.length
                ? `<div class="domain-chip-row">
                    ${entry.individualDomains
                      .map(
                        (domain) =>
                          `<span class="domain-chip">${escapeHtml(domain.label)} ${domain.start}-${domain.end}</span>`
                      )
                      .join("")}
                  </div>`
                : ""
            }
          </div>
        </details>
      `;
      detailPanelEl.innerHTML = `
        <div class="detail-header">
          <div class="detail-header-top">
            <span class="detail-header-name">${escapeHtml(entry.id)}</span>
            <span class="metric-chip">length: ${entry.proteinSequence.length} aa</span>
            <span class="metric-chip ${dotToTone(analysis.validation.lengthStatus)}">CDS: ${escapeHtml(analysis.validation.lengthStatus)}</span>
            <span class="metric-chip">range: ${escapeHtml(formatRange(activeRange))}</span>
            ${(() => {
              const v = analysis.validation;
              if (v.translationStatus === "Match") {
                return `<span class="metric-chip ${dotToTone("Match")}">translation: Match</span>`;
              }
              const deltaPart = v.lengthDelta !== 0 ? ` (CDS transl. ${v.cdsPeptideLength} aa, delta ${v.lengthDelta > 0 ? "+" : ""}${v.lengthDelta})` : "";
              const mismatchPart = v.firstMismatchAa !== null ? `, first mismatch aa ${v.firstMismatchAa}` : "";
              return `<span class="metric-chip ${dotToTone("Mismatch")}" title="${escapeHtml("Translation Mismatch" + deltaPart + mismatchPart)}">translation: Mismatch${mismatchPart}</span>`;
            })()}
          </div>
          <div class="detail-evidence-row">
            ${detailEvidenceCards.map((card) => `
              <span class="detail-evidence-pill ${badgeToneClass(card.cls)}" title="${escapeHtml(card.title)}">
                <span class="detail-evidence-pill-label">${escapeHtml(card.label)}</span>
                <span class="detail-evidence-pill-value">${escapeHtml(card.value)}</span>
              </span>
            `).join("")}
          </div>
        </div>

        ${
          hasStructureModel && state.structureDockRight
            ? `<div class="detail-top-layout">
                 <div class="detail-top-column">
                   ${coordinatesSectionMarkup}
                   ${evidenceSectionMarkup}
                 </div>
                 <div class="detail-top-column detail-top-column-right">${structureDockedMarkup}</div>
               </div>`
            : coordinatesSectionMarkup
        }

        ${hasStructureModel && state.structureDockRight ? "" : evidenceSectionMarkup}

        ${
          hasStructureModel && !state.structureDockRight
            ? structureSectionMarkup
            : ""
        }

        ${
          showCdsIssuePanel
            ? `
              <details class="detail-section" open>
                <summary class="detail-section-header detail-section-summary">
                  <div>
                    <div class="detail-section-index">${cdsIssueSectionIndex}. CDS / Translation Issues</div>
                    <div class="detail-section-copy">Explicit mismatch and truncation view for the current protein or selected construct.</div>
                  </div>
                  <div class="section-chip-row">
                    <span class="pill pill-red">needs review</span>
                    <span class="metric-chip">expand / collapse</span>
                  </div>
                </summary>
                ${renderCdsIssuePanel(entry, analysis, activeRange, construct)}
              </details>
            `
            : ""
        }

        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">${constructSummarySectionIndex}. Construct Summary</div>
              <div class="detail-section-copy">Final selected AA/CDS span, QC summary, and export preview.</div>
            </div>
            <div class="section-chip-row">
              <span class="pill ${toneClass(construct.status)}">${escapeHtml(construct.status)}</span>
              <span class="metric-chip">expand / collapse</span>
            </div>
          </summary>

          <div class="construct-summary-body">
            <div class="construct-summary-grid">
              <div class="summary-card">
                <span>AA range</span>
                <strong>${escapeHtml(formatRange(activeRange))}</strong>
              </div>
              <div class="summary-card">
                <span>CDS range</span>
                <strong>${escapeHtml(formatRange(construct.cdsRange))}</strong>
              </div>
              <div class="summary-card">
                <span>CDS QC</span>
                <strong class="${statusClass(analysis.validation.lengthStatus)}">${escapeHtml(analysis.validation.lengthStatus)}</strong>
              </div>
              <div class="summary-card">
                <span>Translation QC</span>
                ${(() => {
                  const v = analysis.validation;
                  if (v.translationStatus === "Match") {
                    return `<strong class="${statusClass("Match")}">Match</strong>`;
                  }
                  const detail = [
                    v.firstMismatchAa !== null ? `first mismatch: aa ${v.firstMismatchAa}` : null,
                    v.lengthDelta !== 0 ? `CDS encodes ${v.cdsPeptideLength} aa (${v.lengthDelta > 0 ? "+" : ""}${v.lengthDelta})` : null
                  ].filter(Boolean).join(" · ");
                  return `<strong class="${statusClass("Mismatch")}">Mismatch</strong>${detail ? `<small style="display:block;font-weight:normal;color:#888;font-size:10px;margin-top:2px">${escapeHtml(detail)}</small>` : ""}`;
                })()}
              </div>
              <div class="summary-card">
                <span>AA length</span>
                <strong>${aaLen}</strong>
              </div>
              <div class="summary-card">
                <span>NT length</span>
                <strong>${cdsLen}</strong>
              </div>
            </div>

            <div class="output-stack">
              <div>
                <div class="controls-section-label">Peptide sequence</div>
                <pre class="output-pre">${escapeHtml(construct.peptide)}</pre>
              </div>
              <div>
                <div class="controls-section-label">CDS sequence</div>
                <pre class="output-pre">${escapeHtml(construct.cds)}</pre>
              </div>
              <div>
                <div class="controls-section-label">Translated construct</div>
                <pre class="output-pre">${escapeHtml(construct.translation)}</pre>
              </div>
              <div>
                <div class="controls-section-label">FASTA preview</div>
                <pre class="fasta-pre">${escapeHtml(fastaPreview || "Construct is not currently export-ready.")}</pre>
              </div>
            </div>
          </div>
        </details>
      `;

      detailPanelEl.querySelectorAll("[data-candidate-range]").forEach((button) => {
        button.addEventListener("click", () => {
          const candidate = analysis.suggestion.candidates[button.getAttribute("data-candidate-range")];
          if (!candidate) {
            return;
          }
          state.manualRanges[entry.id] = clampRange({
            start: candidate.start,
            end: candidate.end
          }, entry.proteinSequence.length);
          render();
        });
      });

      detailPanelEl.querySelector("[data-toggle-structure-dock]")?.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
        state.structureDockRight = !state.structureDockRight;
        render();
      });

      function applyManualRange(start, end) {
        state.manualRanges[entry.id] = clampRange({ start, end }, entry.proteinSequence.length);
        render();
      }

      function svgClientXToAa(clientX) {
        const svg = detailPanelEl.querySelector("#range-track-svg");
        if (!svg) {
          return 1;
        }

        const labelWidth = Number(svg.dataset.labelWidth);
        const svgWidth = Number(svg.dataset.svgWidth);
        const trackWidth = Number(svg.dataset.trackWidth);
        const proteinLength = Number(svg.dataset.proteinLength);
        const rect = svg.getBoundingClientRect();
        const svgX = ((clientX - rect.left) / rect.width) * svgWidth;
        const fraction = clamp((svgX - labelWidth) / trackWidth, 0, 1);
        return Math.round(1 + fraction * (proteinLength - 1));
      }

      function bindTrackDragging() {
        const svg = detailPanelEl.querySelector("#range-track-svg");
        if (!svg) {
          return;
        }

        svg.querySelectorAll("[data-drag-handle]").forEach((handleEl) => {
          handleEl.addEventListener("mousedown", (event) => {
            event.preventDefault();
            const handle = handleEl.getAttribute("data-drag-handle");
            if (!handle) {
              return;
            }

            svg.classList.add("dragging");

            const onMove = (moveEvent) => {
              const aa = svgClientXToAa(moveEvent.clientX);
              const currentRange = clampRange(
                state.manualRanges[entry.id] ?? getActiveRange(analysis),
                entry.proteinSequence.length
              );

              if (handle === "start") {
                state.manualRanges[entry.id] = clampRange({
                  start: Math.min(aa, currentRange.end - 1),
                  end: currentRange.end
                }, entry.proteinSequence.length);
              } else {
                state.manualRanges[entry.id] = clampRange({
                  start: currentRange.start,
                  end: Math.max(aa, currentRange.start + 1)
                }, entry.proteinSequence.length);
              }
              debounceRender();
            };

            const onUp = () => {
              window.removeEventListener("mousemove", onMove);
              window.removeEventListener("mouseup", onUp);
              detailPanelEl.querySelector("#range-track-svg")?.classList.remove("dragging");
              svg.classList.remove("dragging");
            };

            window.addEventListener("mousemove", onMove);
            window.addEventListener("mouseup", onUp);
          });
        });
      }

      detailPanelEl.querySelector("#range-start-input")?.addEventListener("input", (event) => {
        state.manualRanges[entry.id] = clampRange({
          start: Number(event.target.value),
          end: getActiveRange(analysis).end
        }, entry.proteinSequence.length);
        debounceRender();
      });

      detailPanelEl.querySelector("#range-end-input")?.addEventListener("input", (event) => {
        state.manualRanges[entry.id] = clampRange({
          start: getActiveRange(analysis).start,
          end: Number(event.target.value)
        }, entry.proteinSequence.length);
        debounceRender();
      });

      initializeStructureViewer(entry, activeRange, state.structureDockRight);
      bindTrackDragging();
    }

    function render() {
      const focusedInput = captureFocusedTextInput();
      const analyses = buildAnalyses();
      const visibleAnalyses = analyses.filter((analysis) => {
        const haystack = [
          analysis.entry.id,
          analysis.entry.mergedDomains.map((domain) => domain.label).join(" "),
          analysis.entry.individualDomains.map((domain) => domain.label).join(" ")
        ]
          .join(" ")
          .toLowerCase();
        const matchesSearch = haystack.includes(state.search.toLowerCase());
        const group = analysisGroup(analysis);
        const matchesFilter =
          state.filter === "all"
            ? true
            : state.filter === "good"
              ? group === "Good"
              : state.filter === "range_mismatch"
                ? group === "Range mismatch"
                : group === "CDS error";
        return matchesSearch && matchesFilter;
      });

      if (state.sort !== "id") {
        visibleAnalyses.sort((a, b) => {
          const sa = evidenceSortScore(a.entry, state.sort);
          const sb = evidenceSortScore(b.entry, state.sort);
          return sa !== sb ? sa - sb : a.entry.id.localeCompare(b.entry.id);
        });
      }

      if (!visibleAnalyses.some((analysis) => analysis.entry.id === state.selectedId)) {
        state.selectedId = visibleAnalyses[0]?.entry.id ?? analyses[0]?.entry.id ?? "";
      }

      const selectedAnalysis =
        visibleAnalyses.find((analysis) => analysis.entry.id === state.selectedId) ??
        analyses.find((analysis) => analysis.entry.id === state.selectedId) ??
        analyses[0] ??
        null;

      const batchKey = `${state.search}|${state.filter}|${state.sort}|${state.selectedId}`;
      renderToolbar(selectedAnalysis, analyses);
      renderTabs();

      if (state.view === "review") {
        if (batchKey !== _lastBatchKey) {
          renderBatch(visibleAnalyses, state.selectedId);
          _lastBatchKey = batchKey;
        }
        renderDetail(selectedAnalysis);
      } else if (state.view === "ranges") {
        renderRangesTab(analyses);
      } else {
        renderMetadataTab();
      }

      restoreFocusedTextInput(focusedInput);
    }

    render();
  </script>
</body>
</html>
"""
    return template.replace("__PAYLOAD__", payload_json)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a self-contained HTML construct-design report from protein, CDS, and domain inputs."
    )
    parser.add_argument(
        "--dataset",
        help="Path to dataset.json produced by construct-generate. When provided, raw input flags are ignored.",
    )
    parser.add_argument(
        "--custom-ranges",
        dest="custom_ranges",
        help="Optional TSV/TAB file with additional construct ranges to overlay in the report. Defaults to <input-dir>/custom_ranges.tsv when present, otherwise <input-dir>/custom_ranges.tab.",
    )
    parser.add_argument(
        "--metadata",
        help="Optional metadata TSV/CSV file for the report metadata tab. Defaults to <input-dir>/metadata.tsv, <input-dir>/metadata.csv, <input-dir>/TF_list.csv, repo-root metadata.tsv, repo-root metadata.csv, or repo-root TF_list.csv when present.",
    )
    parser.add_argument(
        "--input-dir",
        help="Base directory used to resolve default input files when explicit paths are omitted. Defaults to the bundled examples/ folder.",
    )
    parser.add_argument(
        "--pep",
        help="Protein FASTA path. Defaults to <input-dir>/proteins.fasta.",
    )
    parser.add_argument(
        "--cds",
        help="CDS FASTA path. Defaults to <input-dir>/cds.fasta.",
    )
    parser.add_argument(
        "--domains",
        help="Optional domain BED path in protein coordinates. Defaults to <input-dir>/domains.bed when present.",
    )
    parser.add_argument(
        "--domains-individual",
        dest="domains_individual",
        help="Optional per-domain BED/TSV path. Defaults to <input-dir>/domains.individual.bed when present, otherwise <input-dir>/domains.individual.tab.",
    )
    parser.add_argument(
        "--evidence-dir",
        dest="evidence_dir",
        help="Optional evidence directory. Defaults to <input-dir>/evidence when present.",
    )
    parser.add_argument(
        "--slop",
        type=int,
        default=DEFAULT_PARAMS["slop"],
        help="Range-calling slop in amino acids.",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=DEFAULT_PARAMS["offset"],
        help="Structure proximity offset in amino acids.",
    )
    parser.add_argument(
        "--min-structured-run",
        dest="min_structured_run",
        type=int,
        default=DEFAULT_PARAMS["minStructuredRun"],
        help="Minimum structured run length used for structure-aware expansion.",
    )
    parser.add_argument(
        "--n-terminal-snap-threshold",
        dest="n_terminal_snap_threshold",
        type=int,
        default=DEFAULT_PARAMS["nTerminalSnapThreshold"],
        help="Snap range ends to the nearest protein terminus when they fall within this threshold.",
    )
    parser.add_argument(
        "--output",
        help="Path to the output HTML file. Defaults to ./report.html in the current working directory.",
    )
    return parser.parse_args(argv)


def parse_generate_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute construct ranges and CDS sequences; write TSV, FASTA, and dataset JSON."
    )
    parser.add_argument("--input-dir")
    parser.add_argument("--pep")
    parser.add_argument("--cds")
    parser.add_argument("--domains")
    parser.add_argument("--domains-individual", dest="domains_individual")
    parser.add_argument("--metadata")
    parser.add_argument("--evidence-dir", dest="evidence_dir")
    parser.add_argument("--slop", type=int, default=DEFAULT_PARAMS["slop"])
    parser.add_argument("--offset", type=int, default=DEFAULT_PARAMS["offset"])
    parser.add_argument(
        "--min-structured-run",
        dest="min_structured_run",
        type=int,
        default=DEFAULT_PARAMS["minStructuredRun"],
    )
    parser.add_argument(
        "--n-terminal-snap-threshold",
        dest="n_terminal_snap_threshold",
        type=int,
        default=DEFAULT_PARAMS["nTerminalSnapThreshold"],
        help="Snap range ends to the nearest protein terminus when they fall within this threshold.",
    )
    parser.add_argument("--output-dir", dest="output_dir", default=".")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    cwd = Path.cwd()
    output_path = (
        resolve_path(cwd, args.output)
        if args.output
        else (cwd / "report.html").resolve()
    )

    if args.dataset:
        dataset_path = resolve_path(cwd, args.dataset)
        if not dataset_path or not dataset_path.exists():
            raise FileNotFoundError(f"Dataset JSON not found: {dataset_path}")
        payload = json.loads(dataset_path.read_text(encoding="utf-8"))
        custom_ranges_path = resolve_path(cwd, args.custom_ranges) if args.custom_ranges else None
        metadata_path = resolve_path(cwd, args.metadata) if args.metadata else None
        if custom_ranges_path and not custom_ranges_path.exists():
            raise FileNotFoundError(f"Custom ranges file not found: {custom_ranges_path}")
        if metadata_path and not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
        payload = append_custom_ranges_to_payload(payload, custom_ranges_path, cwd)
        payload = append_metadata_to_payload(payload, metadata_path, cwd)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        report_from_bundle(payload, output_path)
        print(f"Wrote {output_path}")
        return

    project_root = Path(__file__).resolve().parents[2]
    paths = resolve_cli_inputs(
        cwd=cwd,
        project_root=project_root,
        input_dir_arg=args.input_dir,
        pep_arg=args.pep,
        cds_arg=args.cds,
        domains_arg=args.domains,
        domains_individual_arg=args.domains_individual,
        custom_ranges_arg=args.custom_ranges,
        metadata_arg=args.metadata,
        evidence_dir_arg=args.evidence_dir,
    )

    pep_path = paths["pep_path"]
    cds_path = paths["cds_path"]
    domains_path = paths["domains_path"]
    for required_path in (pep_path, cds_path):
        if not required_path.exists():
            raise FileNotFoundError(f"Required input file not found: {required_path}")
    if domains_path and not domains_path.exists():
        raise FileNotFoundError(f"Optional domain file was requested but not found: {domains_path}")
    for line in format_resolved_input_lines(cwd=cwd, paths=paths):
        print(line)

    params = build_params(
        slop=args.slop,
        offset=args.offset,
        min_structured_run=args.min_structured_run,
        n_terminal_snap_threshold=args.n_terminal_snap_threshold,
    )
    dataset_bundle = load_dataset_bundle(
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        domains_individual_path=paths["domains_individual_path"],
        evidence_root=paths["evidence_root"],
    )
    constructs = compute_constructs(dataset_bundle["entries"], params)
    input_summary = build_input_summary(
        cwd=cwd,
        input_dir=paths["input_dir"],
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        explicit_sources=paths["explicit_sources"],
    )
    payload = build_dataset_payload(
        bundle=dataset_bundle,
        constructs=constructs,
        params=params,
        input_summary=input_summary,
        generated_at=datetime.now(timezone.utc).isoformat(),
    )
    custom_ranges_path = paths["custom_ranges_path"]
    metadata_path = paths["metadata_path"]
    if custom_ranges_path and not custom_ranges_path.exists():
        raise FileNotFoundError(f"Custom ranges file not found: {custom_ranges_path}")
    if metadata_path and not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    payload = append_custom_ranges_to_payload(payload, custom_ranges_path, cwd)
    payload = append_metadata_to_payload(payload, metadata_path, cwd)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_from_bundle(payload, output_path)
    print(f"Wrote {output_path}")
    for line in format_run_summary_lines(dataset_bundle["summary"]):
        print(line)


def generate_main(argv: list[str] | None = None) -> None:
    args = parse_generate_args(argv)
    cwd = Path.cwd()
    project_root = Path(__file__).resolve().parents[2]
    paths = resolve_cli_inputs(
        cwd=cwd,
        project_root=project_root,
        input_dir_arg=args.input_dir,
        pep_arg=args.pep,
        cds_arg=args.cds,
        domains_arg=args.domains,
        domains_individual_arg=args.domains_individual,
        metadata_arg=args.metadata,
        evidence_dir_arg=args.evidence_dir,
    )
    output_dir = resolve_path(cwd, args.output_dir) or cwd
    output_dir.mkdir(parents=True, exist_ok=True)

    pep_path = paths["pep_path"]
    cds_path = paths["cds_path"]
    domains_path = paths["domains_path"]
    for required_path in (pep_path, cds_path):
        if not required_path.exists():
            raise FileNotFoundError(f"Required input file not found: {required_path}")
    if domains_path and not domains_path.exists():
        raise FileNotFoundError(f"Optional domain file was requested but not found: {domains_path}")
    for line in format_resolved_input_lines(cwd=cwd, paths=paths):
        print(line)

    params = build_params(
        slop=args.slop,
        offset=args.offset,
        min_structured_run=args.min_structured_run,
        n_terminal_snap_threshold=args.n_terminal_snap_threshold,
    )
    dataset_bundle = load_dataset_bundle(
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        domains_individual_path=paths["domains_individual_path"],
        evidence_root=paths["evidence_root"],
    )
    constructs = compute_constructs(dataset_bundle["entries"], params)
    input_summary = build_input_summary(
        cwd=cwd,
        input_dir=paths["input_dir"],
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        explicit_sources=paths["explicit_sources"],
    )

    write_constructs_tsv(constructs, output_dir / "constructs.tsv")
    write_constructs_fasta(constructs, output_dir / "constructs.fasta")
    dataset_payload = build_dataset_payload(
        bundle=dataset_bundle,
        constructs=constructs,
        params=params,
        input_summary=input_summary,
    )
    dataset_payload = append_custom_ranges_to_payload(
        dataset_payload,
        paths["custom_ranges_path"],
        cwd,
    )
    dataset_payload = append_metadata_to_payload(
        dataset_payload,
        paths["metadata_path"],
        cwd,
    )
    (output_dir / "dataset.json").write_text(
        json.dumps(dataset_payload, ensure_ascii=False),
        encoding="utf-8",
    )

    for line in format_run_summary_lines(dataset_bundle["summary"]):
        print(line)
    print(f"Wrote {output_dir / 'constructs.tsv'}")
    print(f"Wrote {output_dir / 'constructs.fasta'}")
    print(f"Wrote {output_dir / 'dataset.json'}")


if __name__ == "__main__":
    main()
