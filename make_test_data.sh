set -euo pipefail

N=30
FULL=examples/fulldata
TEST=examples/test

FOLDING=~/ant/gzolotarov/projects/2021_TFevol/folding/results_nvec
CONSERVATION=~/Documents/projects/2025_nvec_motif/probe_design/results/conservation
CONSERVATION_FULL=~/Documents/projects/2025_nvec_motif/probe_design/results/conservation_full
STRUCTURE_DSSP=~/Documents/projects/2025_nvec_motif/probe_design/results/structure

mkdir -p "$TEST/evidence/structures" \
         "$TEST/evidence/conservation" \
         "$TEST/evidence/conservation_full" \
         "$TEST/evidence/structure_dssp"

# all proteins with domain hits
cut -f1 "$FULL/domains.individual.bed" | sort -u > /tmp/domain_ids.txt

# prefer proteins that have a PDB on the folding mount
ls "$FOLDING"/*.alphafold.pdb 2>/dev/null \
  | sed 's/.*\///; s/\.1\.alphafold\.pdb//' \
  | grep -Fx -f /tmp/domain_ids.txt \
  > /tmp/pdb_ids.txt

# take up to N from PDB-covered proteins, pad from the rest if needed
head -n "$N" /tmp/pdb_ids.txt > /tmp/test.ids
HAVE=$(wc -l < /tmp/test.ids)
if [ "$HAVE" -lt "$N" ]; then
    grep -Fxv -f /tmp/test.ids /tmp/domain_ids.txt | head -n $(( N - HAVE )) >> /tmp/test.ids
fi

# sequence subsets
xargs samtools faidx "$FULL/proteins.fasta" < /tmp/test.ids > "$TEST/proteins.fasta"
xargs samtools faidx "$FULL/cds.fasta"      < /tmp/test.ids > "$TEST/cds.fasta"
samtools faidx "$TEST/proteins.fasta"
samtools faidx "$TEST/cds.fasta"

# domain annotations
grep -Fw -f /tmp/test.ids "$FULL/domains.individual.bed" > "$TEST/domains.individual.bed"

# evidence — copy whatever is available for each protein
while IFS= read -r ID; do
    PDB="$FOLDING/${ID}.1.alphafold.pdb"
    TSV="$FOLDING/${ID}.1_plddt_mqc.tsv"
    [ -f "$PDB" ] && cp "$PDB" "$TEST/evidence/structures/${ID}.pdb"
    [ -f "$TSV" ] && cp "$TSV" "$TEST/evidence/structures/${ID}.tsv"

    [ -f "$CONSERVATION/${ID}.out" ]      && cp "$CONSERVATION/${ID}.out"      "$TEST/evidence/conservation/"
    [ -f "$CONSERVATION_FULL/${ID}.out" ] && cp "$CONSERVATION_FULL/${ID}.out"  "$TEST/evidence/conservation_full/"
    [ -f "$STRUCTURE_DSSP/${ID}.out" ]    && cp "$STRUCTURE_DSSP/${ID}.out"     "$TEST/evidence/structure_dssp/"
done < /tmp/test.ids

# summary
echo "Test dataset written to $TEST"
echo "Proteins : $(grep -c '>' "$TEST/proteins.fasta")"
echo "PDB      : $(ls "$TEST/evidence/structures/"*.pdb  2>/dev/null | wc -l)"
echo "pLDDT    : $(ls "$TEST/evidence/structures/"*.tsv  2>/dev/null | wc -l)"
echo "Cons DBD : $(ls "$TEST/evidence/conservation/"*.out       2>/dev/null | wc -l)"
echo "Cons full: $(ls "$TEST/evidence/conservation_full/"*.out  2>/dev/null | wc -l)"
echo "DSSP     : $(ls "$TEST/evidence/structure_dssp/"*.out     2>/dev/null | wc -l)"
