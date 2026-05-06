#!/usr/bin/env bash
set -euo pipefail

N="${N:-30}"
FULL="${FULL:-examples/fulldata}"
TEST="${TEST:-examples/test}"
IDS_FILE=""

usage() {
  cat <<'EOF'
Usage: bash make_test_data.sh [--ids path/to/ids.txt]

Options:
  --ids FILE   Use the protein IDs listed in FILE (one ID per line) instead of
               automatic selection from the full example bundle.
  -h, --help   Show this help text.

Environment overrides:
  N            Number of proteins to select automatically when --ids is absent.
  FULL         Source bundle directory. Default: examples/fulldata
  TEST         Output test bundle directory. Default: examples/test
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --ids)
      [ "$#" -ge 2 ] || { echo "ERROR: --ids requires a file path." >&2; exit 1; }
      IDS_FILE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

PROTEINS="$FULL/proteins.fasta"
CDS="$FULL/cds.fasta"
DOMAINS="$FULL/domains.individual.bed"
CUSTOM_RANGES="$FULL/custom_ranges.tsv"
METADATA="$FULL/metadata.tsv"

STRUCTURES_DIR="$FULL/evidence/structures"
CONSERVATION_DIR="$FULL/evidence/conservation"
CONSERVATION_FULL_DIR="$FULL/evidence/conservation_full"
IUPRED_DIR="$FULL/evidence/iupred"
STRUCTURE_DSSP_DIR="$FULL/evidence/structure_dssp"
STRUCTURE_UNIPROT_DIR="$FULL/evidence/structure_uniprot"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p \
  "$TEST/evidence/structures" \
  "$TEST/evidence/conservation" \
  "$TEST/evidence/conservation_full" \
  "$TEST/evidence/iupred" \
  "$TEST/evidence/structure_dssp" \
  "$TEST/evidence/structure_uniprot"

find "$TEST" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
mkdir -p \
  "$TEST/evidence/structures" \
  "$TEST/evidence/conservation" \
  "$TEST/evidence/conservation_full" \
  "$TEST/evidence/iupred" \
  "$TEST/evidence/structure_dssp" \
  "$TEST/evidence/structure_uniprot"

# all proteins with domain hits
cut -f1 "$DOMAINS" | sort -u > "$TMPDIR/domain_ids.txt"

if [ -n "$IDS_FILE" ]; then
  if [ ! -f "$IDS_FILE" ]; then
    echo "ERROR: --ids file not found: $IDS_FILE" >&2
    exit 1
  fi
  grep -Fvx '' "$IDS_FILE" | sort -u > "$TMPDIR/requested.ids"
  grep -Fx -f "$TMPDIR/requested.ids" "$TMPDIR/domain_ids.txt" > "$TMPDIR/test.ids" || true
  if [ ! -s "$TMPDIR/test.ids" ]; then
    echo "ERROR: none of the requested IDs are present in $DOMAINS" >&2
    exit 1
  fi
else
  # prefer proteins that already have a local structure model in the full example bundle
  find "$STRUCTURES_DIR" -maxdepth 1 -name '*.pdb' -type f \
    | sed 's#^.*/##; s/\.pdb$//' \
    | grep -Fx -f "$TMPDIR/domain_ids.txt" \
    > "$TMPDIR/pdb_ids.txt" || true

  # take up to N from structure-covered proteins, pad from the rest if needed
  head -n "$N" "$TMPDIR/pdb_ids.txt" > "$TMPDIR/test.ids"
  have="$(wc -l < "$TMPDIR/test.ids")"
  if [ "$have" -lt "$N" ]; then
    grep -Fxv -f "$TMPDIR/test.ids" "$TMPDIR/domain_ids.txt" | head -n $(( N - have )) >> "$TMPDIR/test.ids" || true
  fi
fi

# sequence subsets
xargs samtools faidx "$PROTEINS" < "$TMPDIR/test.ids" > "$TEST/proteins.fasta"
xargs samtools faidx "$CDS" < "$TMPDIR/test.ids" > "$TEST/cds.fasta"
samtools faidx "$TEST/proteins.fasta"
samtools faidx "$TEST/cds.fasta"

# domain annotations and per-protein metadata
grep -Fw -f "$TMPDIR/test.ids" "$DOMAINS" > "$TEST/domains.individual.bed"
awk 'BEGIN { FS = OFS = "\t" } NR == FNR { keep[$1] = 1; next } FNR == 1 || ($1 in keep)' \
  "$TMPDIR/test.ids" "$CUSTOM_RANGES" > "$TEST/custom_ranges.tsv"
awk 'BEGIN { FS = OFS = "\t" } NR == FNR { keep[$1] = 1; next } FNR == 1 || ($1 in keep)' \
  "$TMPDIR/test.ids" "$METADATA" > "$TEST/metadata.tsv"

# evidence — copy whatever is available for each protein
while IFS= read -r id; do
  [ -f "$STRUCTURES_DIR/${id}.pdb" ] && cp "$STRUCTURES_DIR/${id}.pdb" "$TEST/evidence/structures/"
  [ -f "$STRUCTURES_DIR/${id}.tsv" ] && cp "$STRUCTURES_DIR/${id}.tsv" "$TEST/evidence/structures/"

  [ -f "$CONSERVATION_DIR/${id}.out" ] && cp "$CONSERVATION_DIR/${id}.out" "$TEST/evidence/conservation/"
  [ -f "$CONSERVATION_FULL_DIR/${id}.out" ] && cp "$CONSERVATION_FULL_DIR/${id}.out" "$TEST/evidence/conservation_full/"
  [ -f "$IUPRED_DIR/${id}.out" ] && cp "$IUPRED_DIR/${id}.out" "$TEST/evidence/iupred/"
  [ -f "$STRUCTURE_DSSP_DIR/${id}.out" ] && cp "$STRUCTURE_DSSP_DIR/${id}.out" "$TEST/evidence/structure_dssp/"
  [ -f "$STRUCTURE_UNIPROT_DIR/${id}.out" ] && cp "$STRUCTURE_UNIPROT_DIR/${id}.out" "$TEST/evidence/structure_uniprot/"
done < "$TMPDIR/test.ids"

count_lines_minus_header() {
  local file="$1"
  if [ ! -f "$file" ]; then
    echo 0
    return
  fi
  awk 'END { print (NR > 0 ? NR - 1 : 0) }' "$file"
}

count_matches() {
  local pattern="$1"
  local dir="${pattern%/*}"
  local name="${pattern##*/}"
  find "$dir" -maxdepth 1 -name "$name" -type f | wc -l
}

echo "Test dataset written to $TEST"
echo "Proteins      : $(grep -c '>' "$TEST/proteins.fasta")"
echo "Custom ranges : $(count_lines_minus_header "$TEST/custom_ranges.tsv")"
echo "Metadata rows : $(count_lines_minus_header "$TEST/metadata.tsv")"
echo "PDB           : $(count_matches "$TEST/evidence/structures/*.pdb")"
echo "pLDDT         : $(count_matches "$TEST/evidence/structures/*.tsv")"
echo "Cons DBD      : $(count_matches "$TEST/evidence/conservation/*.out")"
echo "Cons full     : $(count_matches "$TEST/evidence/conservation_full/*.out")"
echo "IUPred        : $(count_matches "$TEST/evidence/iupred/*.out")"
echo "DSSP          : $(count_matches "$TEST/evidence/structure_dssp/*.out")"
echo "UniProt SS    : $(count_matches "$TEST/evidence/structure_uniprot/*.out")"
if [ -n "$IDS_FILE" ]; then
  echo "ID source     : $IDS_FILE"
else
  echo "ID source     : automatic"
fi
