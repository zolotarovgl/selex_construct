from __future__ import annotations

from construct_report.ranges import (
    build_suggestions,
    extract_structured_runs,
    expand_range_from_nearby_structure,
)


PARAMS = {
    "slop": 10,
    "offset": 15,
    "minStructuredRun": 6,
    "nTerminalSnapThreshold": 15,
}


def _entry(
    seq: str = "A" * 100,
    individual: list[dict[str, int | str]] | None = None,
    merged: list[dict[str, int | str]] | None = None,
    dssp_values: list[str | None] | None = None,
) -> dict[str, object]:
    evidence: dict[str, object] = {}
    if dssp_values is not None:
        evidence["structureDssp"] = {
            "compatible": True,
            "values": dssp_values,
            "coverage": 1.0,
        }
    return {
        "proteinSequence": seq,
        "individualDomains": individual or [],
        "mergedDomains": merged or [],
        "evidence": evidence,
    }


def test_r1_is_domain_envelope() -> None:
    entry = _entry(individual=[{"start": 10, "end": 50, "label": "Foo"}])
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r1"] == {"start": 10, "end": 50}


def test_r2_adds_slop() -> None:
    entry = _entry(individual=[{"start": 20, "end": 60, "label": "Foo"}])
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r2"] == {"start": 10, "end": 70}


def test_r2_clamped_to_protein() -> None:
    entry = _entry(seq="A" * 50, individual=[{"start": 1, "end": 50, "label": "Foo"}])
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r2"] == {"start": 1, "end": 50}


def test_extract_structured_runs_filters_short_segments() -> None:
    track = {"values": [None, "H", "H", "H", None, "E", "E", "E", "E", "E", "E"]}
    runs = extract_structured_runs(track, min_run_length=4)
    assert runs == [{"start": 6, "end": 11}]


def test_r3_expands_into_nearby_structure() -> None:
    values: list[str | None] = [None] * 100
    for index in range(60, 75):
        values[index] = "H"
    entry = _entry(
        individual=[{"start": 20, "end": 60, "label": "Foo"}],
        dssp_values=values,
    )
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r3"]["end"] > suggestions["r2"]["end"]


def test_expand_range_from_nearby_structure_preserves_unrelated_runs() -> None:
    expanded = expand_range_from_nearby_structure(
        {"start": 20, "end": 40},
        structured_runs=[{"start": 60, "end": 70}],
        offset=5,
        protein_length=100,
    )
    assert expanded == {"start": 20, "end": 40}


def test_n_terminal_snap() -> None:
    entry = _entry(seq="A" * 200, individual=[{"start": 20, "end": 80, "label": "Foo"}])
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r3"]["start"] == 1


def test_c_terminal_snap() -> None:
    entry = _entry(seq="A" * 200, individual=[{"start": 120, "end": 190, "label": "Foo"}])
    suggestions = build_suggestions(entry, PARAMS)
    assert suggestions["r3"]["end"] == 200


def test_no_domains_returns_none() -> None:
    suggestions = build_suggestions(_entry(), PARAMS)
    assert suggestions == {"r1": None, "r2": None, "r3": None, "recommended": None}
