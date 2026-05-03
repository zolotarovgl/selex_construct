from __future__ import annotations

import json

from construct_report.pipeline import (
    build_dataset_payload,
    compute_constructs,
    translate_cds,
    write_constructs_fasta,
    write_constructs_tsv,
    write_dataset_json,
)


PARAMS = {
    "slop": 10,
    "offset": 5,
    "minStructuredRun": 4,
    "nTerminalSnapThreshold": 5,
}


def _entry(
    protein_id: str = "P1",
    seq: str = "MAKT" * 25,
    domains: list[dict[str, int | str]] | None = None,
) -> dict[str, object]:
    codon_map = {"M": "ATG", "A": "GCT", "K": "AAA", "T": "ACT"}
    cds = "".join(codon_map.get(residue, "GCT") for residue in seq)
    return {
        "id": protein_id,
        "proteinSequence": seq,
        "cdsSequence": cds,
        "individualDomains": domains if domains is not None else [{"start": 10, "end": 70, "label": "Foo"}],
        "mergedDomains": [],
        "evidence": {},
    }


def test_compute_constructs_returns_one_per_entry_with_domains() -> None:
    constructs = compute_constructs([_entry()], PARAMS)
    assert len(constructs) == 1
    construct = constructs[0]
    assert construct["id"] == "P1"
    assert construct["aaStart"] <= construct["aaEnd"]
    assert construct["cdsLen"] == construct["aaLen"] * 3


def test_compute_constructs_skips_entries_without_domains() -> None:
    constructs = compute_constructs([_entry(domains=[])], PARAMS)
    assert constructs == []


def test_translate_basic() -> None:
    assert translate_cds("ATGGCT") == "MA"


def test_translate_stops_at_stop_codon() -> None:
    assert translate_cds("ATGTAAATGGCT") == "M"


def test_write_tsv(tmp_path) -> None:
    constructs = compute_constructs([_entry()], PARAMS)
    output_path = tmp_path / "constructs.tsv"
    write_constructs_tsv(constructs, output_path)
    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("id\t")
    assert lines[1].startswith("P1\t")


def test_write_fasta_only_includes_ok_constructs(tmp_path) -> None:
    constructs = [
        {
            "id": "P1",
            "recommended": "r3",
            "aaStart": 1,
            "aaEnd": 10,
            "cdsStart": 1,
            "cdsEnd": 30,
            "cds": "ATG" * 10,
            "status": "OK",
        },
        {
            "id": "P2",
            "recommended": "r2",
            "aaStart": 1,
            "aaEnd": 5,
            "cdsStart": 1,
            "cdsEnd": 15,
            "cds": "ATG" * 5,
            "status": "Translation mismatch",
        },
    ]
    output_path = tmp_path / "constructs.fasta"
    write_constructs_fasta(constructs, output_path)
    text = output_path.read_text(encoding="utf-8")
    assert ">P1" in text
    assert ">P2" not in text


def test_build_dataset_payload_shape() -> None:
    bundle = {"entries": [_entry()], "summary": {"keptProteins": 1}}
    constructs = compute_constructs(bundle["entries"], PARAMS)
    payload = build_dataset_payload(bundle, constructs, PARAMS, "test input", generated_at="2026-05-03T00:00:00+00:00")
    assert payload["schemaVersion"] == 1
    assert payload["params"] == PARAMS
    assert len(payload["dataset"]) == 1
    assert len(payload["constructs"]) == 1


def test_write_dataset_json(tmp_path) -> None:
    bundle = {"entries": [_entry()], "summary": {"keptProteins": 1}}
    constructs = compute_constructs(bundle["entries"], PARAMS)
    output_path = tmp_path / "dataset.json"
    write_dataset_json(bundle, constructs, PARAMS, "test input", output_path)
    data = json.loads(output_path.read_text(encoding="utf-8"))
    assert data["schemaVersion"] == 1
    assert data["params"] == PARAMS
    assert "dataset" in data
    assert "constructs" in data
