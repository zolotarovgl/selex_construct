from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


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


def summarize_dataset_for_display(summary: dict[str, Any]) -> str:
    parts = [
        f"pep {summary.get('pepRecords', 0)}",
        f"cds {summary.get('cdsRecords', 0)}",
        f"shared {summary.get('sharedPepCds', 0)}",
        f"kept {summary.get('keptProteins', 0)}",
    ]
    if summary.get("evidenceProteinsAny") or summary.get("evidenceTrackFiles") or summary.get("structureModelFiles"):
        parts.append(f"evidence {summary.get('evidenceProteinsAny', 0)}")
    if summary.get("keptWithStructureModels"):
        parts.append(f"structures {summary.get('keptWithStructureModels', 0)}")
    return " | ".join(parts)


def _format_generated_at(value: str | None) -> str:
    if not value:
        return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return value
    return parsed.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def normalize_report_payload(payload: dict[str, Any]) -> dict[str, Any]:
    dataset_summary = payload.get("datasetSummary") or {}

    if isinstance(dataset_summary, dict) and "display" in dataset_summary and "defaultParams" in payload:
        return payload

    summary_counts = dataset_summary if isinstance(dataset_summary, dict) else {}
    return {
        "schemaVersion": payload.get("schemaVersion", 1),
        "dataset": payload.get("dataset", []),
        "constructs": payload.get("constructs", []),
        "datasetSummary": {
            "counts": summary_counts,
            "display": summarize_dataset_for_display(summary_counts),
        },
        "defaultParams": payload.get("params") or payload.get("defaultParams") or DEFAULT_PARAMS,
        "codonTable": payload.get("codonTable") or CODON_TABLE,
        "inputSummary": payload.get("inputSummary", ""),
        "generatedAt": _format_generated_at(payload.get("generatedAt")),
        "customRanges": payload.get("customRanges", {}),
    }


def render_html(payload_json: str) -> str:
    from construct_report.cli import render_html as legacy_render_html

    return legacy_render_html(payload_json)


def render_payload(payload: dict[str, Any]) -> str:
    normalized = normalize_report_payload(payload)
    payload_json = json.dumps(normalized, separators=(",", ":")).replace("</", "<\\/")
    return render_html(payload_json)


def report_from_bundle(payload: dict[str, Any], output_path: Path) -> None:
    output_path.write_text(render_payload(payload), encoding="utf-8")
