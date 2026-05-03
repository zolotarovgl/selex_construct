from __future__ import annotations

from typing import Any


STRUCTURED_CODES: frozenset[str] = frozenset("H G I E B T".split())


def merge_domains_across_protein(domains: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not domains:
        return []

    sorted_domains = sorted(
        domains,
        key=lambda domain: (domain["start"], domain["end"], str(domain.get("label", ""))),
    )
    labels: list[str] = []
    for domain in sorted_domains:
        label = str(domain.get("label", "")).strip()
        if label and label not in labels:
            labels.append(label)

    return [
        {
            "start": min(domain["start"] for domain in sorted_domains),
            "end": max(domain["end"] for domain in sorted_domains),
            "label": ",".join(labels) or "merged",
        }
    ]


def effective_merged_domains(entry: dict[str, Any]) -> list[dict[str, Any]]:
    individual_domains = entry.get("individualDomains") or []
    if individual_domains:
        return merge_domains_across_protein(individual_domains)

    merged_domains = entry.get("mergedDomains") or []
    if merged_domains:
        return merge_domains_across_protein(merged_domains)

    return []


def envelope_range(domains: list[dict[str, Any]]) -> dict[str, int] | None:
    if not domains:
        return None

    return {
        "start": min(domain["start"] for domain in domains),
        "end": max(domain["end"] for domain in domains),
    }


def normalize_range(range_dict: dict[str, int], protein_length: int) -> dict[str, int]:
    start = max(1, min(int(range_dict["start"]), protein_length))
    end = max(start, min(int(range_dict["end"]), protein_length))
    return {"start": start, "end": end}


def expand_domain_blocks(
    domains: list[dict[str, Any]],
    slop: int,
    protein_length: int,
) -> list[dict[str, Any]]:
    if not domains:
        return []

    expanded_domains = [
        {
            **domain,
            "start": max(1, int(domain["start"]) - slop),
            "end": min(protein_length, int(domain["end"]) + slop),
        }
        for domain in domains
    ]
    return merge_domains_across_protein(expanded_domains)


def preferred_structure_track(entry: dict[str, Any]) -> dict[str, Any] | None:
    evidence = entry.get("evidence") or {}
    if (evidence.get("structureDssp") or {}).get("compatible"):
        return evidence["structureDssp"]
    if (evidence.get("structureUniprot") or {}).get("compatible"):
        return evidence["structureUniprot"]
    return None


def extract_structured_runs(
    track: dict[str, Any] | None,
    min_run_length: int,
) -> list[dict[str, int]]:
    if not track:
        return []

    runs: list[dict[str, int]] = []
    current_start: int | None = None
    values = track.get("values") or []

    for index, value in enumerate(values):
        is_structured = bool(value) and str(value) in STRUCTURED_CODES
        if is_structured and current_start is None:
            current_start = index + 1
        if not is_structured and current_start is not None:
            run = {"start": current_start, "end": index}
            if run["end"] - run["start"] + 1 >= min_run_length:
                runs.append(run)
            current_start = None

    if current_start is not None:
        run = {"start": current_start, "end": len(values)}
        if run["end"] - run["start"] + 1 >= min_run_length:
            runs.append(run)

    return runs


def expand_range_from_nearby_structure(
    base_range: dict[str, int],
    structured_runs: list[dict[str, int]],
    offset: int,
    protein_length: int,
) -> dict[str, int]:
    next_range = dict(base_range)

    for run in structured_runs:
        touches_left = (
            run["start"] <= next_range["start"] + offset
            and run["end"] >= next_range["start"] - offset
        )
        touches_right = (
            run["start"] <= next_range["end"] + offset
            and run["end"] >= next_range["end"] - offset
        )

        if touches_left:
            next_range["start"] = min(next_range["start"], max(1, run["start"] - 1))
        if touches_right:
            next_range["end"] = max(next_range["end"], min(protein_length, run["end"] + 1))

    return normalize_range(next_range, protein_length)


def build_suggestions(entry: dict[str, Any], params: dict[str, Any]) -> dict[str, Any]:
    protein_sequence = entry.get("proteinSequence") or ""
    protein_length = len(protein_sequence)
    r1_segments = effective_merged_domains(entry)
    domain_range = envelope_range(r1_segments)

    if not protein_length or not domain_range:
        return {"r1": None, "r2": None, "r3": None, "recommended": None}

    slop = max(1, int(params["slop"]))
    offset = max(1, int(params["offset"]))
    min_structured_run = max(1, int(params["minStructuredRun"]))
    snap_threshold = max(1, int(params["nTerminalSnapThreshold"]))

    r1 = normalize_range(domain_range, protein_length)
    r2_segments = expand_domain_blocks(r1_segments, slop, protein_length)
    r2 = normalize_range(envelope_range(r2_segments) or domain_range, protein_length)

    structure_track = preferred_structure_track(entry)
    structured_runs = (
        extract_structured_runs(structure_track, min_structured_run)
        if structure_track and structure_track.get("compatible")
        else []
    )
    r3 = expand_range_from_nearby_structure(r2, structured_runs, offset, protein_length)

    if r3["start"] > 1 and r3["start"] < snap_threshold:
        r3 = {**r3, "start": 1}
    if r3["end"] < protein_length and protein_length - r3["end"] < snap_threshold:
        r3 = {**r3, "end": protein_length}

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "recommended": "r3",
    }
