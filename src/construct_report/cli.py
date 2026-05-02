#!/usr/bin/env python3
from __future__ import annotations

import argparse
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
        protein_id, start, end, label, _pfam_id, score = line.split("\t")
        domains.setdefault(protein_id, []).append(
            {
                "start": int(start),
                "end": int(end),
                "label": label,
                "source": "domains.individual.tab",
                "score": float(score),
            }
        )
    return domains


def parse_cases(text: str) -> dict[str, dict[str, str]]:
    lines = clean_lines(text)
    header = lines[0].split("\t")
    rows: dict[str, dict[str, str]] = {}

    for line in lines[1:]:
        values = line.split("\t")
        row = {key: values[index] if index < len(values) else "" for index, key in enumerate(header)}
        rows[row["gene"]] = row

    return rows


def locate_sequence(target: str, query: str) -> int | None:
    if target == query:
        return 1
    match_index = target.find(query)
    return match_index + 1 if match_index >= 0 else None


def parse_conservation_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    if len(lines) < 2:
        raise ValueError("Conservation file must contain at least 2 lines.")

    sequence = lines[0]
    if len(lines) == 2:
        return {"sequence": sequence, "values": [float(value) for value in lines[1].split(",")]}

    alignment = list(lines[1])
    raw_scores = [float(value) for value in lines[2].split(",")]
    values = [score for score, residue in zip(raw_scores, alignment) if residue != "-"]
    return {"sequence": sequence, "values": values}


def parse_structure_file(text: str) -> dict[str, Any]:
    lines = clean_lines(text)
    sequence_index = next((i for i, line in enumerate(lines) if line.endswith("@sequence")), -1)
    structure_index = next((i for i, line in enumerate(lines) if line.endswith("@secondary")), -1)

    if sequence_index < 0 or structure_index < 0:
        raise ValueError("Structure file is missing sequence or secondary blocks.")

    sequence = lines[sequence_index + 1] if sequence_index + 1 < len(lines) else ""
    secondary = lines[structure_index + 1] if structure_index + 1 < len(lines) else ""
    return {"sequence": sequence, "values": list(secondary)}


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

    return {
        "type": "numeric",
        "label": label,
        "sequence": raw_track["sequence"],
        "offset": offset,
        "values": values,
        "coverage": len(raw_track["values"]) / len(protein_sequence),
        "compatible": True,
    }


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


def build_evidence_for_protein(
    protein_id: str,
    protein_sequence: str,
    evidence_files: list[Path],
    structure_files: list[Path],
) -> dict[str, Any]:
    bundle: dict[str, Any] = {}

    # Evidence support is intentionally explicit for now: only these known
    # subdirectory names are loaded into report tracks.
    for path in evidence_files:
        if not path.name.startswith(protein_id):
            continue

        content = path.read_text(encoding="utf-8")
        kind = path.parent.name

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
        elif kind == "structure_uniprot":
            bundle["structureUniprot"] = map_structure_track(
                protein_sequence,
                "UniProt structure",
                parse_structure_file(content),
            )
        elif kind == "structure_dssp":
            bundle["structureDssp"] = map_structure_track(
                protein_sequence,
                "Local DSSP structure",
                parse_structure_file(content),
            )

    matching_structures = sorted(
        (path for path in structure_files if path.name.startswith(protein_id)),
        key=lambda path: (path.stem != protein_id, path.name),
    )
    if matching_structures:
        structure_path = matching_structures[0]
        bundle["structureModel"] = {
            "format": structure_path.suffix.lstrip(".") or "pdb",
            "source": structure_path.name,
            "text": structure_path.read_text(encoding="utf-8"),
        }

    return bundle


def load_dataset(
    pep_path: Path,
    cds_path: Path,
    domains_path: Path | None,
    domains_individual_path: Path | None = None,
    cases_path: Path | None = None,
    evidence_root: Path | None = None,
) -> list[dict[str, Any]]:
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
    cases = (
        parse_cases(cases_path.read_text(encoding="utf-8"))
        if cases_path and cases_path.exists()
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

    entries: list[dict[str, Any]] = []
    for protein_id in sorted(proteins):
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
                ),
                "reference": cases.get(protein_id),
            }
        )

    return entries


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
      --font-ui: "Courier New", monospace;
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
      font-size: 13px;
      line-height: 1.4;
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
      font-size: 14px;
    }

    .app-toolbar-title span {
      font-size: 11px;
      color: var(--muted);
    }

    .app-toolbar-actions {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
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
      grid-template-columns: minmax(260px, 300px) minmax(0, 1fr);
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
      padding: 6px 8px;
      border-bottom: 1px solid var(--border-light);
      background: #f4f4f4;
    }

    .batch-panel-header h2 {
      margin: 0;
      font-size: 12px;
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
      display: block;
      width: 100%;
      padding: 6px 8px;
      border: none;
      border-bottom: 1px solid var(--border-inner);
      border-left: 3px solid transparent;
      background: var(--panel);
      text-align: left;
      cursor: pointer;
      min-height: 48px;
    }

    .review-row:hover {
      background: #f0f4f7;
    }

    .review-row-active {
      border-left-color: var(--accent);
      background: #e8eef4;
    }

    .review-row-top {
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 6px;
      margin-bottom: 2px;
    }

    .review-row-top strong {
      font-size: 12px;
    }

    .review-row-id {
      font-size: 10px;
      color: var(--muted);
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      max-width: 180px;
    }

    .review-row-dots {
      display: flex;
      gap: 4px;
      align-items: center;
    }

    .status-dot {
      display: inline-block;
      width: 8px;
      height: 8px;
      border-radius: 50%;
    }

    .dot-green { background: var(--good); }
    .dot-red { background: var(--bad); }
    .dot-amber { background: var(--amber); }
    .dot-gray { background: #aaa; }

    .detail-column {
      padding: 10px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      min-width: 0;
    }

    .detail-header {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      align-items: center;
      padding-bottom: 2px;
      border-bottom: 1px solid var(--border-light);
    }

    .detail-header-name {
      font-size: 14px;
      font-weight: 700;
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
    }

    .coordinates-extra-wrap {
      padding: 0 10px 10px;
      display: flex;
      flex-direction: column;
      gap: 8px;
      background: #fff;
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

    .structure-panel {
      border: 1px solid var(--border-light);
      background: #fff;
    }

    .structure-panel-header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 8px;
      padding: 8px 10px;
      border-bottom: 1px solid var(--border-light);
      background: #fafafa;
    }

    .structure-panel-meta {
      display: flex;
      gap: 4px;
      flex-wrap: wrap;
      justify-content: flex-end;
    }

    .structure-viewer {
      height: 340px;
      background: #fff;
    }

    .structure-fallback {
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 340px;
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

      .controls-row {
        grid-template-columns: 1fr;
      }

      .app-toolbar {
        flex-direction: column;
        align-items: stretch;
      }
    }
  </style>
</head>
<body>
  <main class="app-shell">
    <div id="toolbar" class="app-toolbar"></div>
    <div class="workspace">
      <aside id="batch-panel" class="panel batch-panel"></aside>
      <section id="detail-panel" class="panel detail-column"></section>
    </div>
  </main>

  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <script id="report-data" type="application/json">__PAYLOAD__</script>
  <script>
    const payload = JSON.parse(document.getElementById("report-data").textContent);
    const dataset = payload.dataset;
    const defaultParams = payload.defaultParams;
    const codonTable = payload.codonTable;
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
    const structuredCodes = new Set(["H", "G", "I", "E", "B", "T"]);
    const initialEntry =
      dataset.find((entry) => entry.evidence?.structureModel?.text) ??
      dataset.find((entry) => Object.keys(entry.evidence || {}).length > 0) ??
      dataset[0] ??
      null;

    const state = {
      params: { ...defaultParams },
      search: "",
      filter: "all",
      selectedId: initialEntry?.id ?? "",
      manualRanges: {},
      showCoordinateDetails: false
    };

    const toolbarEl = document.getElementById("toolbar");
    const batchPanelEl = document.getElementById("batch-panel");
    const detailPanelEl = document.getElementById("detail-panel");

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

      return {
        lengthStatus,
        translationStatus,
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
      return dataset.map((entry) => ({
        entry,
        suggestion: buildSuggestions(entry, state.params),
        validation: validateProtein(entry),
        referenceRange: parseRangeString(entry.reference?.picked_range ?? "")
      }));
    }

    function getActiveRange(analysis) {
      const recommendedRange =
        (analysis.suggestion.recommendedKey
          ? analysis.suggestion.candidates[analysis.suggestion.recommendedKey]
          : null) ?? analysis.referenceRange;

      return clampRange(
        state.manualRanges[analysis.entry.id] ??
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

    const BROWSER_LABEL_W = 124;
    const BROWSER_SVG_W = 1040;
    const BROWSER_TRACK_W = BROWSER_SVG_W - BROWSER_LABEL_W;

    function renderTrackViewer(entry, candidates, activeRange) {
      const LABEL_W = BROWSER_LABEL_W;
      const SVG_W = BROWSER_SVG_W;
      const TRACK_W = BROWSER_TRACK_W;
      const GAP = 3;
      const mergedDomains = effectiveMergedDomains(entry);
      const hasIndividualDomains = entry.individualDomains.length > 0;
      const individualLaneHeight = 10;
      const individualLaneGap = 4;

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
            8 + // top padding
            8 + // individual row label
            4 + // gap after individual label
            individualLaneCount * individualLaneHeight +
            Math.max(0, individualLaneCount - 1) * individualLaneGap +
            8 + // gap before divider / merged section
            8 + // merged row label
            4 + // gap after merged label
            12 + // merged row height
            8 + // gap before contribution lane
            10 + // contribution lane
            6 // bottom padding
          )
        : 54;
      const rows = [
        { label: "", height: 30, type: "ruler" },
        { label: "Domains", height: domainRowHeight, type: "domains" }
      ];

      function xPos(pos) {
        if (entry.proteinSequence.length <= 1) {
          return LABEL_W;
        }
        return LABEL_W + ((pos - 1) / (entry.proteinSequence.length - 1)) * TRACK_W;
      }

      const candidateByKey = Object.fromEntries(candidates.map((candidate) => [candidate.key, candidate]));
      const r1 = candidateByKey.r1 ?? null;
      const r2 = candidateByKey.r2 ?? null;
      const r3 = candidateByKey.r3 ?? null;
      const r1Segments = r1?.segments?.length ? r1.segments : [];
      const r2Segments = r2?.segments?.length ? r2.segments : [];
      const tickInterval =
        entry.proteinSequence.length > 400 ? 100 : entry.proteinSequence.length > 150 ? 50 : 25;

      const ticks = [1];
      for (let pos = tickInterval; pos < entry.proteinSequence.length; pos += tickInterval) {
        ticks.push(pos);
      }
      if (ticks[ticks.length - 1] !== entry.proteinSequence.length) {
        ticks.push(entry.proteinSequence.length);
      }

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

      function pushLabeledDomainBox(kind, domain, x, y, width, height, color, opacity, fontSize) {
        const clipId = `${kind}-domain-clip-${domainClipCounter++}`;
        svg.push(
          `<g><title>${escapeHtml(domain.label)} ${domain.start}-${domain.end}</title>` +
            `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="${color}" opacity="${opacity}" rx="1.5" stroke="#000" stroke-width="1" />` +
            `<defs><clipPath id="${clipId}"><rect x="${x + 1}" y="${y + 1}" width="${Math.max(0, width - 2)}" height="${Math.max(0, height - 2)}" rx="1.5" /></clipPath></defs>` +
            `<text x="${x + 3}" y="${y + height / 2 + fontSize / 3 - 0.5}" font-size="${fontSize}" fill="#111" clip-path="url(#${clipId})">${escapeHtml(domain.label)}</text>` +
          `</g>`
        );
      }

      const domainRowTop = trackTops[1];
      const individualLabelY = domainRowTop + 14;
      const individualTrackTop = domainRowTop + 20;
      const individualBlockHeight =
        individualLaneCount * individualLaneHeight +
        Math.max(0, individualLaneCount - 1) * individualLaneGap;
      const mergedLabelY = hasIndividualDomains ? individualTrackTop + individualBlockHeight + 16 : domainRowTop + 14;
      const mergedY = hasIndividualDomains ? mergedLabelY + 4 : domainRowTop + 7;
      const mergedHeight = hasIndividualDomains ? 12 : 18;

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
          row.type === "domains" && hasIndividualDomains
            ? `<text x="${LABEL_W - 4}" y="${mergedLabelY}" text-anchor="end" font-size="7" fill="#666">merged</text>`
            : "",
          index > 0 ? `<line x1="0" y1="${y}" x2="${SVG_W}" y2="${y}" stroke="#d0d0d0" stroke-width="1" />` : "",
          `<line x1="${LABEL_W}" y1="${y}" x2="${LABEL_W}" y2="${y + row.height}" stroke="#aaa" stroke-width="1" />`
        );
      });

      svg.push(
        `<rect x="${selectedStartX}" y="0" width="${Math.max(0, selectedEndX - selectedStartX)}"` +
          ` height="${totalHeight}" fill="rgba(44,62,80,0.10)" />`
      );

      const rulerY = trackTops[0];
      const rulerBaseline = rulerY + rows[0].height - 2;
      svg.push(`<line x1="${LABEL_W}" y1="${rulerBaseline}" x2="${SVG_W}" y2="${rulerBaseline}" stroke="#bbb" stroke-width="1" />`);
      ticks.forEach((pos) => {
        const x = xPos(pos);
        svg.push(
          `<line x1="${x}" y1="${rulerBaseline - 5}" x2="${x}" y2="${rulerBaseline}" stroke="#999" stroke-width="1" />`,
          `<text x="${x}" y="${rulerBaseline - 7}" text-anchor="middle" font-size="8" fill="#666">${pos}</text>`
        );
      });

      candidates.forEach((candidate) => {
        const color = candidateColors[candidate.key] ?? "#888";
        svg.push(
          `<text x="${xPos(candidate.start) + 3}" y="${rulerY + 10}" font-size="8" fill="${color}" font-weight="700">` +
            `${candidate.key}</text>`
        );
      });

      const contributionY = mergedY + mergedHeight + 7;
      const contributionHeight = 10;
      if (!mergedDomains.length && !entry.individualDomains.length) {
        svg.push(
          `<text x="${LABEL_W + 10}" y="${mergedY + 14}" font-size="10" fill="#888">no domain BED ranges</text>`
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
        svg.push(
          `<line x1="${LABEL_W}" y1="${individualTrackTop + individualBlockHeight + 8}" x2="${SVG_W}" y2="${individualTrackTop + individualBlockHeight + 8}" stroke="#d4d4d4" stroke-width="1" />`
        );
      }

      mergedDomains.forEach((domain, index) => {
        const x = xPos(domain.start);
        const width = Math.max(4, xPos(domain.end) - x);
        const color = trackPalette[index % trackPalette.length];
        pushLabeledDomainBox("merged", domain, x, mergedY, width, mergedHeight, color, 0.72, hasIndividualDomains ? 7 : 8);
      });

      if (r2Segments.length) {
        r2Segments.forEach((segment) => {
          const segmentX = xPos(segment.start);
          const segmentWidth = Math.max(3, xPos(segment.end) - segmentX);
          svg.push(
            `<rect x="${segmentX}" y="${contributionY}" width="${segmentWidth}" height="${contributionHeight}"` +
              ` fill="rgba(230,126,34,0.35)" stroke="#e67e22" stroke-width="1" rx="1" />`
          );
        });
      }

      if (r1Segments.length) {
        r1Segments.forEach((segment) => {
          const segmentX = xPos(segment.start);
          const segmentWidth = Math.max(3, xPos(segment.end) - segmentX);
          svg.push(
            `<rect x="${segmentX}" y="${contributionY}" width="${segmentWidth}" height="${contributionHeight}"` +
              ` fill="rgba(231,76,60,0.20)" stroke="#e74c3c" stroke-width="1" rx="1" />`
          );
        });
      }

      if (r2 && r3) {
        const r3Segments = [
          r3.start < r2.start ? { start: r3.start, end: r2.start } : null,
          r3.end > r2.end ? { start: r2.end, end: r3.end } : null
        ].filter(Boolean);

        r3Segments.forEach((segment) => {
          const segmentX = xPos(segment.start);
          const segmentWidth = Math.max(3, xPos(segment.end) - segmentX);
          svg.push(
            `<rect x="${segmentX}" y="${contributionY}" width="${segmentWidth}" height="${contributionHeight}"` +
              ` fill="rgba(39,174,96,0.35)" stroke="#27ae60" stroke-width="1" rx="1" />`
          );
          if (segmentWidth > 34) {
            svg.push(
              `<text x="${segmentX + 4}" y="${contributionY + 8}" font-size="7.5" fill="#1f7a45" font-weight="700">r3 ext</text>`
            );
          }
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

      svg.push(
        `<line x1="${selectedStartX}" y1="0" x2="${selectedStartX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="3" />`,
        `<line x1="${selectedEndX}" y1="0" x2="${selectedEndX}" y2="${totalHeight}" stroke="${selectedRangeColor}" stroke-width="3" />`,
        `<circle cx="${selectedStartX}" cy="10" r="4" fill="${selectedRangeColor}" />`,
        `<circle cx="${selectedEndX}" cy="10" r="4" fill="${selectedRangeColor}" />`
      );

      svg.push(
        `<line x1="${selectedStartX}" y1="0" x2="${selectedStartX}" y2="${totalHeight}"` +
          ` stroke="transparent" stroke-width="18" data-drag-handle="start" />`,
        `<line x1="${selectedEndX}" y1="0" x2="${selectedEndX}" y2="${totalHeight}"` +
          ` stroke="transparent" stroke-width="18" data-drag-handle="end" />`
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
          height: 22
        });
      }
      if (entry.evidence.structureUniprot?.compatible) {
        availableTracks.push({
          type: "structure",
          label: "UniProt SS",
          coverage: entry.evidence.structureUniprot.coverage,
          values: entry.evidence.structureUniprot.values,
          height: 22
        });
      }
      if (entry.evidence.conservationFull?.compatible) {
        availableTracks.push({
          type: "numeric",
          label: "Cons full",
          coverage: entry.evidence.conservationFull.coverage,
          values: entry.evidence.conservationFull.values,
          height: 40,
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
          height: 40,
          fill: "rgba(39,174,96,0.28)",
          stroke: "#1f7a45"
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

      const rows = [{ type: "ruler", label: "", height: 24 }, ...availableTracks];

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
            ? `<text x="${LABEL_W - 6}" y="${y + row.height / 2 - 1}" text-anchor="end" font-size="9" fill="#444" font-weight="700">${escapeHtml(row.label)}</text>`
            : "",
          row.type !== "ruler"
            ? `<text x="${LABEL_W - 6}" y="${y + row.height / 2 + 8}" text-anchor="end" font-size="8" fill="#777">${Math.round(row.coverage * 100)}%</text>`
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
            svg.push(
              `<rect x="${x}" y="${y + 4}" width="${width}" height="${track.height - 8}" fill="${color}" rx="1" />`
            );
          });
        } else {
          const trackTop = y + 4;
          const trackHeight = track.height - 8;
          const maxValue = Math.max(
            1,
            ...track.values.filter((value) => value !== null && value !== undefined).map(Number)
          );
          const areaPath = buildAreaPath(track.values, xPos, trackTop, trackHeight, maxValue);
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
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-numeric"></span>numeric conservation track</span>
            <span class="browser-legend-item"><span class="browser-legend-swatch browser-legend-range"></span>selected construct span</span>
          </div>
        </div>
      `;
    }

    function renderStructurePanel(entry, activeRange) {
      const model = entry.evidence.structureModel;
      if (!model?.text) {
        return `
          <div class="structure-panel">
            <div class="structure-panel-header">
              <div>
                <div class="controls-section-label">Structure viewer</div>
                <div class="detail-section-copy">Local structure model with the selected construct range highlighted.</div>
              </div>
              <div class="structure-panel-meta">
                <span class="metric-chip">no model</span>
              </div>
            </div>
            <div class="structure-fallback">No matching PDB model was found under the evidence <code>structures/</code> directory for this protein.</div>
          </div>
        `;
      }

      return `
        <div class="structure-panel">
          <div class="structure-panel-header">
            <div>
              <div class="controls-section-label">Structure viewer</div>
              <div class="detail-section-copy">Local structure model with the selected construct range highlighted.</div>
            </div>
            <div class="structure-panel-meta">
              <span class="metric-chip">${escapeHtml(model.source)}</span>
              <span class="metric-chip">AA ${escapeHtml(formatRange(activeRange))}</span>
            </div>
          </div>
          <div
            id="structure-viewer-${sanitizeDomId(entry.id)}"
            class="structure-viewer"
            data-structure-viewer
          ></div>
        </div>
      `;
    }

    function initializeStructureViewer(entry, activeRange) {
      const viewerEl = detailPanelEl.querySelector("[data-structure-viewer]");
      const model = entry.evidence.structureModel;
      if (!viewerEl || !model?.text) {
        return;
      }

      if (!window.$3Dmol || typeof window.$3Dmol.createViewer !== "function") {
        viewerEl.innerHTML = `
          <div class="structure-fallback">
            The 3D viewer library did not load. Open the report with internet access to fetch 3Dmol.js, or vendor the script locally for offline viewing.
          </div>
        `;
        return;
      }

      try {
        const viewer = window.$3Dmol.createViewer(viewerEl, { backgroundColor: "white" });
        viewer.addModel(model.text, model.format || "pdb");
        viewer.setStyle({}, { cartoon: { color: "#d0d0d0" } });
        viewer.setStyle(
          { resi: `${activeRange.start}-${activeRange.end}` },
          {
            cartoon: { color: selectedRangeColor },
            stick: { radius: 0.18, color: "#e67e22" }
          }
        );
        viewer.zoomTo();
        viewer.render();
      } catch (error) {
        const message = error instanceof Error ? error.message : String(error);
        viewerEl.innerHTML = `
          <div class="structure-fallback">
            Failed to render the local PDB model: ${escapeHtml(message)}
          </div>
        `;
      }
    }

    function renderSequencePanel(entry, activeRange) {
      const aaStart = activeRange.start - 1;
      const aaEnd = activeRange.end;
      const cdsStart = aaStart * 3;
      const cdsEnd = aaEnd * 3;

      const aaBefore = escapeHtml(entry.proteinSequence.slice(0, aaStart));
      const aaSelected = escapeHtml(entry.proteinSequence.slice(aaStart, aaEnd));
      const aaAfter = escapeHtml(entry.proteinSequence.slice(aaEnd));

      const cdsBefore = escapeHtml(entry.cdsSequence.slice(0, cdsStart));
      const cdsSelected = escapeHtml(entry.cdsSequence.slice(cdsStart, cdsEnd));
      const cdsAfter = escapeHtml(entry.cdsSequence.slice(cdsEnd));

      return `
        <div class="seq-panel">
          <div class="seq-track-label">AA · ${activeRange.start}-${activeRange.end} (${aaSelected.length} aa)</div>
          <pre class="seq-pre"><span class="seq-ctx">${aaBefore}</span><mark class="seq-mark">${aaSelected}</mark><span class="seq-ctx">${aaAfter}</span></pre>
          <div class="seq-track-label">CDS · ${cdsStart + 1}-${cdsEnd} (${cdsSelected.length} nt)</div>
          <pre class="seq-pre"><span class="seq-ctx">${cdsBefore}</span><mark class="seq-mark">${cdsSelected}</mark><span class="seq-ctx">${cdsAfter}</span></pre>
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
          const recommended =
            analysis.suggestion.recommendedKey &&
            analysis.suggestion.candidates[analysis.suggestion.recommendedKey];

          if (!recommended) {
            return [];
          }

          const itemConstruct = buildConstruct(analysis.entry, recommended);
          if (itemConstruct.status !== "OK") {
            return [];
          }

          return [
            `>${analysis.entry.id}|${recommended.key}|aa=${formatRange(recommended)}\\n${itemConstruct.cds}`
          ];
        })
        .join("\\n");

      toolbarEl.innerHTML = `
        <div class="app-toolbar-title">
          <strong>Construct Design HTML Report</strong>
          <span>Required inputs: protein FASTA, CDS FASTA, optional domain BED</span>
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
        </div>
        <div class="batch-toolbar">
          <input id="search-input" value="${escapeHtml(state.search)}" placeholder="Search sequence ID or domain…" aria-label="Search">
          <select id="filter-input" aria-label="Filter">
            <option value="all"${state.filter === "all" ? " selected" : ""}>All</option>
            <option value="good"${state.filter === "good" ? " selected" : ""}>Good only</option>
            <option value="attention"${state.filter === "attention" ? " selected" : ""}>Attention</option>
          </select>
        </div>
        <div class="review-list">
          ${
            visibleAnalyses.length
              ? visibleAnalyses
                  .map((analysis) => {
                    const reference = analysis.entry.reference ?? {};
                    const dots = [
                      reference.status_structure ?? "NA",
                      reference.status_range ?? "NA",
                      analysis.validation.lengthStatus
                    ];
                    const active = analysis.entry.id === selectedId;
                    const domainSummary =
                      analysis.entry.mergedDomains.map((domain) => domain.label).join(", ") ||
                      "no domain BED ranges";

                    return `
                      <button type="button" class="review-row ${active ? "review-row-active" : ""}" data-select-id="${escapeHtml(analysis.entry.id)}">
                        <div class="review-row-top">
                          <strong>${escapeHtml(analysis.entry.id)}</strong>
                          <span class="review-row-dots">
                            ${dots
                              .map((status) => `<span class="status-dot ${dotClass(status)}" title="${escapeHtml(status)}"></span>`)
                              .join("")}
                          </span>
                        </div>
                        <div class="review-row-id">${escapeHtml(domainSummary)}</div>
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
        render();
      });

      batchPanelEl.querySelector("#filter-input")?.addEventListener("change", (event) => {
        state.filter = event.target.value;
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
      const coordinateDetailsLabel = state.showCoordinateDetails ? "hide lower details" : "show lower details";
      detailPanelEl.innerHTML = `
        <div class="detail-header">
          <span class="detail-header-name">${escapeHtml(entry.id)}</span>
          <span class="metric-chip">length: ${entry.proteinSequence.length} aa</span>
          <span class="metric-chip">CDS: <span class="${statusClass(analysis.validation.lengthStatus)}">&nbsp;${escapeHtml(analysis.validation.lengthStatus)}</span></span>
          <span class="metric-chip">range: ${escapeHtml(formatRange(activeRange))}</span>
          <span class="metric-chip">translation: <span class="${statusClass(analysis.validation.translationStatus)}">&nbsp;${escapeHtml(analysis.validation.translationStatus)}</span></span>
        </div>

        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">1. Coordinates and Sequence</div>
              <div class="detail-section-copy">Protein coordinates with optional range controls and AA/CDS highlights below.</div>
            </div>
            <div class="section-chip-row">
              <span class="metric-chip">selected: ${escapeHtml(formatRange(activeRange))}</span>
              <span class="metric-chip">${aaLen} aa</span>
              <span class="metric-chip">${cdsLen} nt</span>
              <button type="button" class="coordinates-extra-toggle" data-toggle-coordinate-details>${coordinateDetailsLabel}</button>
              <span class="metric-chip">expand / collapse</span>
            </div>
          </summary>

          <div class="track-viewer-wrap">
            ${renderTrackViewer(entry, candidateRanges, activeRange)}
          </div>
          ${
            state.showCoordinateDetails
              ? `
                <div class="coordinates-extra-wrap">
                  <div class="track-legend">
                    <span class="track-legend-item"><span class="track-legend-swatch track-legend-r1"></span>r1 merged domain span</span>
                    <span class="track-legend-item"><span class="track-legend-swatch track-legend-r2"></span>r2 r1 plus slop</span>
                    <span class="track-legend-item"><span class="track-legend-swatch track-legend-r3"></span>r3 extra beyond r2</span>
                  </div>
                  <details class="track-explainer">
                    <summary>What are r1 / r2 / r3?</summary>
                    <div class="track-explainer-body">
                      <p><strong>r1</strong> is built by merging all available domain hits for the protein into one continuous span.</p>
                      <p><strong>r2</strong> is that <strong>r1</strong> span expanded by ±${state.params.slop} aa.</p>
                      <p><strong>r3</strong> starts from <strong>r2</strong> and expands further when a compatible structure track has structured runs of at least ${state.params.minStructuredRun} aa within ${state.params.offset} aa of either edge.</p>
                      <p>If the resulting start lands below residue ${state.params.nTerminalSnapThreshold} but is still greater than 1, the range snaps to residue 1.</p>
                      <p>So in practice, <strong>r3</strong> is the structure-aware version of <strong>r2</strong>, meant to avoid cutting too close to nearby structured elements.</p>
                    </div>
                  </details>
                  <div class="track-drag-note">Drag the black range boundaries directly on the coordinate plot.</div>
                  ${renderSequencePanel(entry, activeRange)}
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
                      <div class="params-panel">
                        <div class="controls-section-label">Suggestion parameters</div>
                        <div class="params-grid">
                          <div class="params-field">
                            <label>Slop
                              <input type="number" data-param-key="slop" value="${state.params.slop}">
                            </label>
                          </div>
                          <div class="params-field">
                            <label>Offset
                              <input type="number" data-param-key="offset" value="${state.params.offset}">
                            </label>
                          </div>
                          <div class="params-field">
                            <label>Min struct run
                              <input type="number" data-param-key="minStructuredRun" value="${state.params.minStructuredRun}">
                            </label>
                          </div>
                          <div class="params-field">
                            <label>N-term snap
                              <input type="number" data-param-key="nTerminalSnapThreshold" value="${state.params.nTerminalSnapThreshold}">
                            </label>
                          </div>
                        </div>
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
                </div>
              `
              : ""
          }
        </details>

        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">2. Evidence Tracks</div>
              <div class="detail-section-copy">Mapped browser tracks and optional local structure models aligned to protein coordinates.</div>
            </div>
            <div class="metric-chip">expand / collapse</div>
          </summary>
          <div class="metadata-body">
            ${renderEvidenceBrowser(entry, activeRange)}
            ${renderStructurePanel(entry, activeRange)}
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

        <details class="detail-section" open>
          <summary class="detail-section-header detail-section-summary">
            <div>
              <div class="detail-section-index">3. Construct Summary</div>
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
                <strong class="${statusClass(analysis.validation.translationStatus)}">${escapeHtml(analysis.validation.translationStatus)}</strong>
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
          state.manualRanges[entry.id] = { start: candidate.start, end: candidate.end };
          render();
        });
      });

      detailPanelEl.querySelector("[data-toggle-coordinate-details]")?.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
        state.showCoordinateDetails = !state.showCoordinateDetails;
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
                applyManualRange(Math.min(aa, currentRange.end - 1), currentRange.end);
              } else {
                applyManualRange(currentRange.start, Math.max(aa, currentRange.start + 1));
              }
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
        const nextStart = Number(event.target.value);
        applyManualRange(nextStart, getActiveRange(analysis).end);
      });

      detailPanelEl.querySelector("#range-end-input")?.addEventListener("input", (event) => {
        const nextEnd = Number(event.target.value);
        applyManualRange(getActiveRange(analysis).start, nextEnd);
      });

      detailPanelEl.querySelectorAll("[data-param-key]").forEach((input) => {
        input.addEventListener("input", (event) => {
          const key = input.getAttribute("data-param-key");
          const value = Math.max(1, Number(event.target.value) || 1);
          state.params[key] = value;
          render();
        });
      });

      initializeStructureViewer(entry, activeRange);
      bindTrackDragging();
    }

    function render() {
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
        const reference = analysis.entry.reference;
        const matchesFilter =
          state.filter === "all"
            ? true
            : state.filter === "good"
              ? reference?.status_range === "Good"
              : reference?.status_range !== "Good";
        return matchesSearch && matchesFilter;
      });

      if (!visibleAnalyses.some((analysis) => analysis.entry.id === state.selectedId)) {
        state.selectedId = visibleAnalyses[0]?.entry.id ?? analyses[0]?.entry.id ?? "";
      }

      const selectedAnalysis =
        visibleAnalyses.find((analysis) => analysis.entry.id === state.selectedId) ??
        analyses.find((analysis) => analysis.entry.id === state.selectedId) ??
        analyses[0] ??
        null;

      renderToolbar(selectedAnalysis, analyses);
      renderBatch(visibleAnalyses, state.selectedId);
      renderDetail(selectedAnalysis);
    }

    render();
  </script>
</body>
</html>
"""
    return template.replace("__PAYLOAD__", payload_json)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a self-contained HTML construct-design report from protein, CDS, and domain inputs."
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
        help="Optional per-domain TSV path. Defaults to <input-dir>/domains.individual.tab when present.",
    )
    parser.add_argument(
        "--cases",
        help="Optional cases TSV path. Defaults to <input-dir>/cases.tsv when present.",
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
        help="Snap starts to residue 1 when they fall below this threshold.",
    )
    parser.add_argument(
        "--output",
        help="Path to the output HTML file. Defaults to ./report.html in the current working directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[2]
    cwd = Path.cwd()
    input_dir = (
        resolve_path(cwd, args.input_dir)
        if args.input_dir
        else (project_root / "examples").resolve()
    )
    primary_explicit = any(value is not None for value in (args.pep, args.cds, args.domains))
    pep_path = resolve_path(cwd, args.pep) or (input_dir / "proteins.fasta")
    cds_path = resolve_path(cwd, args.cds) or (input_dir / "cds.fasta")
    if args.domains:
        domains_path = resolve_path(cwd, args.domains)
    elif primary_explicit:
        domains_path = None
    else:
        default_domains = input_dir / "domains.bed"
        domains_path = default_domains if default_domains.exists() else None
    domains_individual_path = (
        resolve_path(cwd, args.domains_individual)
        if args.domains_individual
        else (None if primary_explicit else input_dir / "domains.individual.tab")
    )
    cases_path = (
        resolve_path(cwd, args.cases)
        if args.cases
        else (None if primary_explicit else input_dir / "cases.tsv")
    )
    evidence_root = (
        resolve_path(cwd, args.evidence_dir)
        if args.evidence_dir
        else (None if primary_explicit else input_dir / "evidence")
    )
    output_path = (
        resolve_path(cwd, args.output)
        if args.output
        else (cwd / "report.html").resolve()
    )

    for required_path in (pep_path, cds_path):
        if not required_path.exists():
            raise FileNotFoundError(f"Required input file not found: {required_path}")
    if domains_path and not domains_path.exists():
        raise FileNotFoundError(f"Optional domain file was requested but not found: {domains_path}")

    params = {
        "slop": max(1, args.slop),
        "offset": max(1, args.offset),
        "minStructuredRun": max(1, args.min_structured_run),
        "nTerminalSnapThreshold": max(1, args.n_terminal_snap_threshold),
    }

    dataset = load_dataset(
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=domains_path,
        domains_individual_path=domains_individual_path,
        cases_path=cases_path,
        evidence_root=evidence_root,
    )
    explicit_sources = any(
        value is not None
        for value in (
            args.pep,
            args.cds,
            args.domains,
            args.domains_individual,
            args.cases,
            args.evidence_dir,
        )
    )
    input_summary = (
        " | ".join(
            [
                f"pep={format_path_for_display(pep_path, cwd)}",
                f"cds={format_path_for_display(cds_path, cwd)}",
                f"domains={format_path_for_display(domains_path, cwd) if domains_path else 'none'}",
            ]
        )
        if explicit_sources
        else format_path_for_display(input_dir, cwd)
    )

    payload = {
        "dataset": dataset,
        "defaultParams": params,
        "codonTable": CODON_TABLE,
        "inputSummary": input_summary,
        "generatedAt": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    }
    payload_json = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    html_text = render_html(payload_json)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html_text, encoding="utf-8")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
