# construct_report — Agent Instructions

This document describes the purpose, inputs, GUI layout, and internal logic of `construct_report` for an LLM agent picking up work on this codebase.

---

## Purpose

`construct_report` generates recommended construct boundaries, CDS exports, and a **self-contained HTML report** for reviewing protein construct boundaries. The typical use case is HT-SELEX construct design: given a set of transcription factor proteins (here *Nematostella vectensis*), the tool suggests which amino-acid range to express recombinantly, integrating domain annotations, sequence conservation, and AlphaFold structure data.

Most orchestration now lives across `src/construct_report/cli.py`, `ranges.py`, `pipeline.py`, and `report.py`; the report step still uses the embedded HTML/JavaScript template renderer. The HTML output remains self-contained and requires no server.

---

## Inputs

### Required

| Flag | Description |
|------|-------------|
| `--pep` | Protein FASTA. IDs must match CDS FASTA exactly. |
| `--cds` | CDS FASTA. One CDS per protein, ID-matched to `--pep`. |

### Optional

| Flag | Description |
|------|-------------|
| `--domains` | Merged domain BED in amino-acid coordinates (`protein_id  start  end  label`). |
| `--domains-individual` | Per-domain-hit BED/TSV (4 or 6 columns: `protein_id  start  end  label [pfam_id  score]`). If provided and non-empty, only proteins appearing in this file are included in the report. |
| `--cases` | TSV with human-curated review annotations per protein (see Cases Table below). |
| `--evidence-dir` | Root directory containing evidence subdirectories (see Evidence Layout below). |
| `--dataset` | Report-only mode for `construct-report`; path to `dataset.json` produced by `construct-generate`. |
| `--output` | Output HTML path. Defaults to `./report.html`. |

### Range-Calling Parameters

| Flag | Default | Effect |
|------|---------|--------|
| `--slop` | 50 | AA added on both sides of r1 to produce r2. |
| `--offset` | 15 | Proximity window (AA) for structure-aware r2→r3 expansion. |
| `--min-structured-run` | 6 | Minimum continuous structured segment length for r3 expansion. |
| `--n-terminal-snap-threshold` | 15 | Snap r3 start to residue 1 when it falls within this many AA of the N-terminus. |

---

## Evidence Directory Layout

```
evidence/
  conservation/          # per-residue conservation scores, domain alignment
  conservation_full/     # per-residue conservation scores, full-length alignment
  iupred/                # local IUPred disorder scores (.out)
  structure_dssp/        # DSSP secondary structure strings from local AlphaFold models
  structure_uniprot/     # UniProt secondary structure strings
  structures/            # AlphaFold PDB models (.pdb) and pLDDT score tables (.tsv)
```

Files are matched to proteins by prefix: a file is used for protein `ID` if `file.name.startswith(ID)`. Malformed or empty files emit a WARNING to stderr and are skipped — they do not crash the run.

The current report code consumes `conservation/`, `conservation_full/`, `iupred/`, `structure_dssp/`, `structure_uniprot/`, and `structures/`.

### Evidence File Formats

**`conservation/` and `conservation_full/`** — two or three lines:
```
SEQUENCE
[ALIGNMENT-WITH-GAPS]    # optional; gaps stripped before score mapping
0.61,0.85,0.79,...       # comma-separated floats
```

#### Detailed provenance for `conservation/` and `conservation_full/`

Both conservation tracks come from the upstream `../../probe_design` analysis, but they are produced from different alignments and represent slightly different biological contexts.

**`conservation/`**

- Source path in the upstream workflow: `../../probe_design/results/conservation/<ID>.out`
- Purpose: per-residue conservation on the older TFevol alignment that was already clipped around the DNA-binding domain or family-defining region
- Typical usage in the report: a domain-centric conservation view that is often useful for seeing how strongly the annotated core region is conserved

The upstream recipe in `../../probe_design/README.md` builds this track by:

1. Building an `ID -> orthology group` mapping from TFevol `gene_trees/*groups.csv`.
2. Looking up the corresponding TFevol alignment in `.../results/align/<HG>.aln.fasta`.
3. Extracting the target sequence from that alignment.
4. Running `../../probe_design/scripts/aln_parse.py --remove-gaps` so that only non-gap positions in the target are retained.

That means `conservation/` is based on an existing alignment that may already have been clipped or filtered for phylogenetic purposes, rather than a fresh full-length realignment.

**`conservation_full/`**

- Source path in the upstream workflow: `../../probe_design/results/conservation_full/<ID>.out`
- Purpose: full-length per-residue conservation over a fresh realignment of all proteins from the same orthology group as the target
- Typical usage in the report: broad construct-design context, especially for evaluating whether suggested termini extend into well-conserved or poorly conserved flanking sequence

This track was introduced because the clipped TFevol alignment is not ideal for judging full-length construct boundaries. The upstream note in `../../probe_design/README.md` explicitly states that a better conservation profile should come from realigning whole sequences from the same orthology group, rather than using the already clipped alignment.

The upstream generation logic for `conservation_full/` is:

1. Concatenate the TFevol `*groups.csv` files to obtain a global `ID -> orthology group` lookup.
2. For a given target `ID`, identify all proteins assigned to the same orthology group.
3. Fetch those full-length protein sequences from the merged Nematostella protein FASTA with `samtools faidx`.
4. Write them to `tmp/realign/<ID>.fasta`.
5. Realign that orthology-group FASTA with `famsa`, producing `tmp/realign/<ID>.aln`.
6. Extract the ungapped target sequence itself into the output file.
7. Run `../../probe_design/scripts/aln_parse.py -i tmp/realign/<ID>.aln -id <ID>` and append the computed scores to the same output file.

The command pattern documented upstream is:

```bash
PROBE=../../probe_design
TMP_DIR=$PROBE/tmp/realign
OUTDIR=$PROBE/results/conservation_full

# Fetch the orthology-group FASTA for each target ID first
famsa -t 12 "$TMP_DIR/<ID>.fasta" "$TMP_DIR/<ID>.aln"
samtools faidx "$TMP_DIR/<ID>.fasta" <ID> | bioawk -c fastx '{print $seq}' > "$OUTDIR/<ID>.out"
python "$PROBE/scripts/aln_parse.py" -i "$TMP_DIR/<ID>.aln" -id <ID> >> "$OUTDIR/<ID>.out"
```

##### Scoring function

`../../probe_design/scripts/aln_parse.py` contains two candidate scoring functions:

- `shannon_entropy(column)`
- `conservation_score(column)`

Despite the argument parser description mentioning Shannon entropy, the active implementation uses `conservation_score(column)`, not entropy. The script currently computes:

```text
most_common_residue_count / column_length
```

for every alignment column, rounded to 2 decimals. So the score is a simple majority-frequency conservation score:

- `1.00` means all aligned residues in that column are identical
- lower values indicate more heterogeneous columns
- higher values therefore correspond to stronger conservation

##### On-disk file layout

`construct_report` accepts both 2-line and 3-line conservation files because both exist in the upstream history.

Two-line form:

```text
SEQUENCE
0.61,0.85,0.79,...
```

Three-line form:

```text
SEQUENCE
ALIGNMENT-WITH-GAPS
0.61,0.85,0.79,...
```

Interpretation:

- line 1: ungapped target protein sequence
- line 2: aligned target sequence from the MSA, including `-` gap characters
- line 3: one conservation value per alignment column

When `construct_report` loads a 3-line conservation file, it removes values at positions where line 2 contains a gap before mapping the scores back to the ungapped target sequence. This mirrors the upstream R logic in `../../probe_design/pick_ranges.R`, where `.load_conservation()` strips scores for `aln == '-'`.

##### Relationship to the report

- `conservation/` is rendered as the `Cons DBD` numeric track
- `conservation_full/` is rendered as the `Cons full` numeric track
- both are optional; missing or empty files are skipped with warnings
- matching is by protein ID prefix, not by an external manifest

For the bundled example datasets, the files under:

- `examples/fulldata/evidence/conservation/`
- `examples/fulldata/evidence/conservation_full/`
- `examples/test/evidence/conservation/`
- `examples/test/evidence/conservation_full/`

were copied from the corresponding upstream `../../probe_design/results/...` directories for proteins present in the example FASTAs.

**`iupred/`** — one tab-separated line per protein:
```
SEQUENCE<TAB>0.2764,0.2680,0.2680,...
```

- Column 1: the ungapped protein sequence passed to IUPred.
- Column 2: comma-separated per-residue IUPred scores in sequence order.
- The number of scores must equal the sequence length.

These files were produced locally from `../../iupred/results/*.out`, based on the command sequence documented in `../../iupred/README.md`. The local wrapper splits the input FASTA into one sequence per file, runs:
```bash
python "$IUPRED" "$FILE" long
```
and then collapses the raw per-residue output into the compact on-disk format used here:
```bash
python "$IUPRED" "$FILE" long \
  | grep -v '#' \
  | awk '{s=s $2; sc=sc $3 ","} END{print s "\t" substr(sc,1,length(sc)-1)}'
```

For the bundled examples, only files whose basename matches a protein ID in the example FASTA were copied:
- `examples/fulldata/evidence/iupred/`: 657 matching proteins
- `examples/test/evidence/iupred/`: 81 matching proteins

**`structure_dssp/` and `structure_uniprot/`** — FASTA-like:
```
>protein_id@sequence
MAKTL...
>protein_id@secondary
---HHHHHTTSSPP---
```

**`structures/<id>.pdb`** — standard PDB format. Used for 3D viewer.

**`structures/<id>.tsv`** — pLDDT score table from AlphaFold multimer output:
```
Positions  rank_0  rank_1  rank_2  rank_3  rank_4
1          47.36   ...
2          46.87   ...
```
Only `rank_0` is loaded. Positions are 1-based model residue indices, mapped to full-protein coordinates using the structure model's `proteinStart` offset.

### Cases Table (`cases.tsv`)

Tab-separated with header. Key columns used by the report:

| Column | Values | Used for |
|--------|--------|----------|
| `gene` | protein ID | row key |
| `status_structure` | Full / Shorter / UniprotBetter / Missing | reference annotation |
| `status_range` | Good / Bad | filter dropdown (Good only / Attention) |
| `picked_range` | e.g. `107-320` | fallback active range when no manual override |
| `picked_range_name` | r1 / r2 / r3 | which candidate was chosen |

---

## Range Logic

Three candidate ranges are computed per protein:

| Range | Definition |
|-------|------------|
| `r1` | Merged span of all domain hits (individual or merged BED) |
| `r2` | r1 expanded by ±`slop` AA, clamped to protein length |
| `r3` | r2 expanded further where structure tracks have continuous structured runs (≥`min-structured-run` AA) within `offset` AA of the r2 boundary. If the resulting start is within `n-terminal-snap-threshold` AA of position 1, it snaps to 1. |

`r3` is the default recommended range. The user can override it by dragging the range handles or editing the numeric inputs in the detail panel.

The **active range** for a protein is:
1. Manual override (`state.manualRanges[id]`), if set by the user.
2. Otherwise, `r3` (or `r2`/`r1` if `r3` is unavailable).
3. Fallback: the `picked_range` from `cases.tsv`.

---

## Two-Step Workflow

`construct-report` still accepts raw inputs directly, but the codebase now also supports a split generate/report workflow.

### Step 1 — Generate constructs

```bash
construct-generate \
  --pep proteins.fasta \
  --cds cds.fasta \
  --domains-individual domains.individual.bed \
  --evidence-dir evidence/ \
  --output-dir out/
```

This writes:

- `out/constructs.tsv` — per-protein ranges, selected CDS spans, and status
- `out/constructs.fasta` — CDS FASTA for all constructs whose translation matches
- `out/dataset.json` — full data bundle for the report step

### Step 2 — Generate HTML report

```bash
construct-report --dataset out/dataset.json --output report.html
```

Or keep the original single-step behaviour:

```bash
construct-report \
  --pep proteins.fasta \
  --cds cds.fasta \
  --output report.html
```

---

## GUI Layout

The report has three panels: **Toolbar** (top), **Batch list** (left), **Detail view** (right).

### Toolbar

- Title and dataset summary (record counts, input paths, generation timestamp).
- **Export selected FASTA**: downloads CDS for the currently active range of the selected protein. Respects manual range edits.
- **Export all FASTA**: downloads CDS for all proteins using each protein's current active range (including any manual overrides).

### Batch List (left panel)

A scrollable list of all proteins in the report. Each row shows:
- **Protein ID** (bold).
- **Evidence badges** (top-right of row): four colored square badges indicating what data is loaded for that protein:
- **Evidence badges** (top-right of row): colored square badges indicating what data is loaded for that protein:

  | Badge | Meaning | Green | Amber | Gray |
  |-------|---------|-------|-------|------|
  | `3D` | AlphaFold PDB model | present | — | absent |
  | `SS` | Secondary structure (DSSP) | DSSP, full coverage | DSSP, partial coverage OR UniProt-only | absent |
  | `CN` | Conservation track | any conservation loaded | — | absent |
  | `ID` | IUPred disorder track | present | — | absent |
  | `pL` | pLDDT scores | present | — | absent |

- **Domain summary** (below ID): comma-separated individual domain labels, or merged domain labels, or "no domains".
- **Left border color**: driven by `status_range` from `cases.tsv` (green = Good, red = Bad, amber = Shorter/UniprotBetter). Active row shows accent blue, overriding the range color.

**Controls above the list:**
- Search box: filters by protein ID or domain label.
- Filter dropdown: All / Good only / Attention (based on `cases.tsv` `status_range`).
- Sort dropdown: ID (alphabetical) / evidence (total badge score) / 3D+SS / conservation / pLDDT. Sorting is descending by evidence richness within the selected category.

### Detail View (right panel)

Shows the full review interface for the selected protein. Sections:

#### Header chips
- Protein length (AA).
- CDS status chip (colored: green = MatchSTOP/Match, red = mismatch, amber = other).
- Active range coordinates.
- Translation status chip.

#### Section 1 — Coordinates and Sequence

- **Range track SVG**: visual ruler with draggable start/end handles for r1, r2, r3, and the active range. Individual and merged domain blocks shown.
- **Range inputs**: numeric start/end fields; debounced 80 ms.
- **Candidate buttons** (r1 / r2 / r3): click to snap the active range to that candidate.
- **Parameters panel** (collapsible): slop, offset, min-structured-run, n-terminal-snap-threshold; changes trigger debounced re-render and invalidate the `buildAnalyses` cache.
- **Sequence display** (collapsible): AA sequence with selected range highlighted; CDS sequence below.

#### Section 2 — Evidence Browser

SVG track browser showing per-residue evidence aligned to the full protein sequence. Track order (top to bottom):
1. Ruler (position ticks).
2. DSSP secondary structure (if compatible).
3. UniProt secondary structure (if compatible).
4. Full-length conservation (blue area, if compatible).
5. Domain conservation / DBD (green area, if compatible).
6. IUPred disorder (gray area, fixed 0–1 scale, if compatible).
7. pLDDT scores (amber area, fixed 0–100 scale, if compatible).

The active range is highlighted as a shaded rectangle with vertical boundary lines.

#### Section 3 — 3D Structure Viewer

Rendered only if a `.pdb` file was loaded for the protein. Uses `3Dmol.js` (loaded from CDN; requires internet). The selected construct range is highlighted in red on the cartoon model.

**Performance note:** the 3Dmol viewer is cached per protein ID. Switching proteins recreates the WebGL context; changing only the range updates colors without recreating the viewer.

#### Section 4 — Construct Output

- Construct status pill (green = OK, red = error).
- CDS length and translation check results.
- FASTA preview of the selected construct (protein ID, AA range, CDS range, wrapped CDS sequence).

---

## Internal Architecture

The main runtime pieces are:

| Component | Role |
|-----------|------|
| `cli.py` | CLI entry points, path resolution, input loading, backward-compatible single-step flow |
| `ranges.py` | Python port of the `r1` / `r2` / `r3` range logic |
| `pipeline.py` | Generate step: construct calling plus TSV / FASTA / JSON writing |
| `report.py` | Report payload normalization and HTML report writing |
| `parse_fasta` | Parse protein / CDS FASTA |
| `parse_merged_domains` / `parse_individual_domains` | Parse domain BED files |
| `parse_conservation_file` | Parse 2- or 3-line conservation format; raises `ValueError` on empty/malformed |
| `parse_structure_file` | Parse DSSP/UniProt SS format |
| `parse_pdb_model` | Extract residue numbers and sequence from PDB |
| `parse_plddt_file` | Parse AlphaFold pLDDT TSV (rank_0 column) |
| `map_numeric_track` / `map_structure_track` / `map_plddt_track` | Map evidence onto full-protein coordinate array |
| `build_structure_model_mapping` | Align PDB model to full protein using DSSP/UniProt offset or sequence search |
| `build_evidence_for_protein` | Assemble all evidence for one protein; skips bad files with stderr warning |
| `load_dataset_bundle` | Load all inputs; return entries + summary dict |
| `render_html` | Embed data as JSON + emit the full HTML/CSS/JS template used by the report step |
| `main` / `generate_main` | Report and generate CLI entry points |

**JavaScript state** (inside the HTML):
- `state.params` — range-calling parameters (slop, offset, etc.)
- `state.search` / `state.filter` / `state.sort` — batch list controls
- `state.selectedId` — currently displayed protein
- `state.manualRanges` — user-edited ranges keyed by protein ID
- `_structureCache` — cached 3Dmol viewer per entry ID (avoids WebGL recreation on re-render)
- `_analysesCache` — memoized `buildAnalyses()` result keyed by params JSON

**Render flow:** `render()` → `buildAnalyses()` (memoized) → filter + sort → `renderToolbar()` + `renderBatch()` (skipped if batch key unchanged) + `renderDetail()`. Range drags and param inputs are debounced (80 ms) before triggering `render()`.

---

## Adding New Evidence Types

New evidence types are **not auto-discovered**. To add one, extend these four locations in `cli.py`:

1. `build_evidence_for_protein()` — add a new `elif kind == "my_type":` branch that parses and maps the file into `bundle["myKey"]`.
2. `load_dataset_bundle()` — scan for new file extensions if needed (currently scans `*.out`, `*.pdb`, `*.tsv`).
3. `renderEvidenceBrowser()` (JS) — add a new `availableTracks.push(...)` block referencing `entry.evidence.myKey`.
4. `evidenceBadges()` (JS) — add or update the badge logic to reflect presence/quality of the new track.
