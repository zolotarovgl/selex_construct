from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from construct_report.ranges import build_suggestions


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


def translate_cds(sequence: str) -> str:
    residues: list[str] = []
    normalized = str(sequence or "").upper()

    for index in range(0, len(normalized) - 2, 3):
        codon = normalized[index : index + 3]
        residue = CODON_TABLE.get(codon, "X")
        if residue == "*":
            break
        residues.append(residue)

    return "".join(residues)


def compute_constructs(
    entries: list[dict[str, Any]],
    params: dict[str, Any],
) -> list[dict[str, Any]]:
    constructs: list[dict[str, Any]] = []

    for entry in entries:
        suggestion = build_suggestions(entry, params)
        recommended_key = suggestion["recommended"]
        if recommended_key is None:
            continue

        active_range = suggestion[recommended_key]
        aa_start = int(active_range["start"])
        aa_end = int(active_range["end"])
        cds_start = (aa_start - 1) * 3 + 1
        cds_end = aa_end * 3

        cds_sequence = entry.get("cdsSequence") or ""
        protein_sequence = entry.get("proteinSequence") or ""
        peptide = protein_sequence[aa_start - 1 : aa_end]
        cds_slice = cds_sequence[cds_start - 1 : min(cds_end, len(cds_sequence))]
        translation = translate_cds(cds_slice)

        if cds_end > len(cds_sequence):
            status = "Range exceeds CDS length"
        elif translation == peptide:
            status = "OK"
        else:
            status = "Translation mismatch"

        constructs.append(
            {
                "id": entry["id"],
                "r1": suggestion["r1"],
                "r2": suggestion["r2"],
                "r3": suggestion["r3"],
                "recommended": recommended_key,
                "aaStart": aa_start,
                "aaEnd": aa_end,
                "cdsStart": cds_start,
                "cdsEnd": cds_end,
                "aaLen": aa_end - aa_start + 1,
                "cdsLen": len(cds_slice),
                "peptide": peptide,
                "cds": cds_slice,
                "translation": translation,
                "status": status,
            }
        )

    return constructs


def _format_range(range_dict: dict[str, Any] | None) -> str:
    if not range_dict:
        return "NA"
    return f"{range_dict['start']}-{range_dict['end']}"


def write_constructs_tsv(constructs: list[dict[str, Any]], path: Path) -> None:
    header = [
        "id",
        "r1",
        "r2",
        "r3",
        "recommended",
        "aa_start",
        "aa_end",
        "cds_start",
        "cds_end",
        "aa_len",
        "cds_len",
        "status",
        "cds",
    ]

    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for construct in constructs:
            row = [
                construct["id"],
                _format_range(construct["r1"]),
                _format_range(construct["r2"]),
                _format_range(construct["r3"]),
                construct["recommended"],
                str(construct["aaStart"]),
                str(construct["aaEnd"]),
                str(construct["cdsStart"]),
                str(construct["cdsEnd"]),
                str(construct["aaLen"]),
                str(construct["cdsLen"]),
                construct["status"],
                construct["cds"],
            ]
            handle.write("\t".join(row) + "\n")


def write_constructs_fasta(constructs: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for construct in constructs:
            if construct["status"] != "OK":
                continue

            header = (
                f">{construct['id']}|{construct['recommended']}"
                f"|aa={construct['aaStart']}-{construct['aaEnd']}"
                f"|cds={construct['cdsStart']}-{construct['cdsEnd']}"
            )
            handle.write(header + "\n")

            cds = construct["cds"]
            for index in range(0, len(cds), 60):
                handle.write(cds[index : index + 60] + "\n")


def build_dataset_payload(
    bundle: dict[str, Any],
    constructs: list[dict[str, Any]],
    params: dict[str, Any],
    input_summary: str,
    generated_at: str | None = None,
) -> dict[str, Any]:
    return {
        "schemaVersion": 1,
        "generatedAt": generated_at or datetime.now(timezone.utc).isoformat(),
        "params": params,
        "inputSummary": input_summary,
        "dataset": bundle["entries"],
        "datasetSummary": bundle["summary"],
        "constructs": constructs,
    }


def write_dataset_json(
    bundle: dict[str, Any],
    constructs: list[dict[str, Any]],
    params: dict[str, Any],
    input_summary: str,
    path: Path,
) -> None:
    payload = build_dataset_payload(
        bundle=bundle,
        constructs=constructs,
        params=params,
        input_summary=input_summary,
    )
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
