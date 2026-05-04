from __future__ import annotations

import json

from construct_report.cli import (
    generate_main,
    load_dataset_bundle,
    parse_custom_ranges,
    parse_iupred_file,
    parse_metadata_table,
)


def test_parse_custom_ranges_supports_constructs_tsv_header() -> None:
    text = "\n".join(
        [
            "id\trecommended\taa_start\taa_end",
            "P1\tr3\t10\t120",
            "P2\tr2\t5\t90",
        ]
    )
    ranges = parse_custom_ranges(text)
    assert ranges["P1"] == [{"start": 10, "end": 120, "label": "r3"}]
    assert ranges["P2"] == [{"start": 5, "end": 90, "label": "r2"}]


def test_parse_custom_ranges_supports_range_column() -> None:
    text = "\n".join(
        [
            "gene\tpicked_range\tpicked_range_name",
            "P1\t20-80\tlegacy",
            "P1\t25-75\tmanual",
        ]
    )
    ranges = parse_custom_ranges(text)
    assert ranges["P1"] == [
        {"start": 20, "end": 80, "label": "legacy"},
        {"start": 25, "end": 75, "label": "manual"},
    ]


def test_parse_custom_ranges_skips_na_picked_ranges() -> None:
    text = "\n".join(
        [
            "gene\tpicked_range\tpicked_range_name",
            "P1\t20-80\tr3",
            "P2\tNA\tr2",
            "P3\t\tmanual",
        ]
    )
    ranges = parse_custom_ranges(text)
    assert ranges == {
        "P1": [{"start": 20, "end": 80, "label": "r3"}],
    }


def test_generate_main_picks_up_input_dir_custom_ranges(tmp_path, capsys) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    evidence_dir = input_dir / "evidence" / "iupred"
    evidence_dir.mkdir(parents=True)

    protein_seq = "M" + ("A" * 29)
    cds_seq = "ATG" + ("GCT" * 29)

    (input_dir / "proteins.fasta").write_text(f">P1\n{protein_seq}\n", encoding="utf-8")
    (input_dir / "cds.fasta").write_text(f">P1\n{cds_seq}\n", encoding="utf-8")
    (input_dir / "domains.individual.bed").write_text(
        "P1\t5\t20\tHomeobox\n",
        encoding="utf-8",
    )
    (input_dir / "custom_ranges.tsv").write_text(
        "gene\tpicked_range\tpicked_range_name\nP1\t4-22\tlegacy\n",
        encoding="utf-8",
    )
    (evidence_dir / "P1.out").write_text(
        f"{protein_seq}\t" + ",".join(["0.10"] * len(protein_seq)) + "\n",
        encoding="utf-8",
    )

    generate_main(["--input-dir", str(input_dir), "--output-dir", str(output_dir)])

    payload = json.loads((output_dir / "dataset.json").read_text(encoding="utf-8"))
    assert payload["customRanges"] == {
        "P1": [{"start": 4, "end": 22, "label": "legacy"}],
    }
    stdout = capsys.readouterr().out
    assert "Picked up from input directory:" in stdout
    assert "custom ranges:" in stdout
    assert "evidence found:" in stdout
    assert "iupred: 1 .out file" in stdout


def test_parse_iupred_file_supports_compact_local_format() -> None:
    parsed = parse_iupred_file("MAKT\t0.10,0.20,0.30,0.40\n")
    assert parsed == {
        "sequence": "MAKT",
        "values": [0.10, 0.20, 0.30, 0.40],
        "maxValue": 1.0,
    }


def test_parse_metadata_table_moves_gene_column_first_and_filters() -> None:
    text = "\n".join(
        [
            "rank\t\tgene\tgene_name\tzscore",
            "1\tDnvec_1\tP1\tAlpha\t12.5",
            "2\tDnvec_2\tP2\tBeta\t8.0",
        ]
    )
    table = parse_metadata_table(text, allowed_ids={"P2"})

    assert [column["key"] for column in table["columns"][:3]] == ["gene", "rank", "column_2"]
    assert table["idKey"] == "gene"
    assert table["rows"] == [
        {
            "gene": "P2",
            "rank": "2",
            "column_2": "Dnvec_2",
            "gene_name": "Beta",
            "zscore": "8.0",
        }
    ]


def test_load_dataset_bundle_loads_iupred_track(tmp_path) -> None:
    pep_path = tmp_path / "proteins.fasta"
    cds_path = tmp_path / "cds.fasta"
    evidence_root = tmp_path / "evidence"
    iupred_dir = evidence_root / "iupred"
    iupred_dir.mkdir(parents=True)

    pep_path.write_text(">P1\nMAKT\n", encoding="utf-8")
    cds_path.write_text(">P1\nATGGCTAAAACT\n", encoding="utf-8")
    (iupred_dir / "P1.out").write_text("MAKT\t0.10,0.20,0.30,0.40\n", encoding="utf-8")

    bundle = load_dataset_bundle(
        pep_path=pep_path,
        cds_path=cds_path,
        domains_path=None,
        evidence_root=evidence_root,
    )

    entry = bundle["entries"][0]
    iupred = entry["evidence"]["iupred"]
    assert iupred["compatible"] is True
    assert iupred["label"] == "IUPred disorder"
    assert iupred["values"] == [0.10, 0.20, 0.30, 0.40]
    assert iupred["coverage"] == 1.0
    assert iupred["maxValue"] == 1.0


def test_generate_main_picks_up_input_dir_metadata(tmp_path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()

    protein_seq = "M" + ("A" * 29)
    cds_seq = "ATG" + ("GCT" * 29)

    (input_dir / "proteins.fasta").write_text(f">P1\n{protein_seq}\n", encoding="utf-8")
    (input_dir / "cds.fasta").write_text(f">P1\n{cds_seq}\n", encoding="utf-8")
    (input_dir / "domains.individual.bed").write_text(
        "P1\t5\t20\tHomeobox\n",
        encoding="utf-8",
    )
    (input_dir / "TF_list.csv").write_text(
        "\n".join(
            [
                "rank\t\tgene\tgene_name\tzscore",
                "1\tDnvec_1\tP1\tAlpha\t12.5",
                "2\tDnvec_2\tP2\tBeta\t8.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    generate_main(["--input-dir", str(input_dir), "--output-dir", str(output_dir)])

    payload = json.loads((output_dir / "dataset.json").read_text(encoding="utf-8"))
    assert payload["metadataTable"]["idKey"] == "gene"
    assert [column["key"] for column in payload["metadataTable"]["columns"][:3]] == [
        "gene",
        "rank",
        "column_2",
    ]
    assert payload["metadataTable"]["rows"] == [
        {
            "gene": "P1",
            "rank": "1",
            "column_2": "Dnvec_1",
            "gene_name": "Alpha",
            "zscore": "12.5",
        }
    ]
