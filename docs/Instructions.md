# construct_report — Agent Instructions

This document describes the purpose, inputs, GUI layout, and internal logic of `construct_report` for an LLM agent picking up work on this codebase.

---

## Purpose

`construct_report` generates a **self-contained HTML report** for reviewing protein construct boundaries. The typical use case is HT-SELEX construct design: given a set of transcription factor proteins (here *Nematostella vectensis*), the tool suggests which amino-acid range to express recombinantly, integrating domain annotations, sequence conservation, and AlphaFold structure data.

All logic — data parsing, range calculation, and the interactive GUI — lives in a single file: `src/construct_report/cli.py`. The HTML output embeds all data and JavaScript inline, requiring no server.

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
  structure_dssp/        # DSSP secondary structure strings from local AlphaFold models
  structure_uniprot/     # UniProt secondary structure strings
  structures/            # AlphaFold PDB models (.pdb) and pLDDT score tables (.tsv)
```

Files are matched to proteins by prefix: a file is used for protein `ID` if `file.name.startswith(ID)`. Malformed or empty files emit a WARNING to stderr and are skipped — they do not crash the run.

### Evidence File Formats

**`conservation/` and `conservation_full/`** — two or three lines:
```
SEQUENCE
[ALIGNMENT-WITH-GAPS]    # optional; gaps stripped before score mapping
0.61,0.85,0.79,...       # comma-separated floats
```

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

  | Badge | Meaning | Green | Amber | Gray |
  |-------|---------|-------|-------|------|
  | `3D` | AlphaFold PDB model | present | — | absent |
  | `SS` | Secondary structure (DSSP) | DSSP, full coverage | DSSP, partial coverage OR UniProt-only | absent |
  | `CN` | Conservation track | any conservation loaded | — | absent |
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
6. pLDDT scores (amber area, fixed 0–100 scale, if compatible).

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

All code lives in `src/construct_report/cli.py`:

| Function | Role |
|----------|------|
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
| `render_html` | Embed data as JSON + emit the full HTML/CSS/JS template |
| `main` | CLI entry point; path resolution, validation, orchestration |

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
