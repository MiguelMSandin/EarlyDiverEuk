#!/bin/bash

# First, remove identical sequences
# We will take the longest sequences as database
ID="0.97"
THREADS="8"
OUT="identicals.tsv"
LIST="identicals.list"

vsearch --threads $THREADS --usearch_global "pr2_otus99.fasta" --db "pacbio.fasta" --id "$ID" --blast6out "$OUT"

# Now we select those PR2 OTUs that have a 100% similarity to the PacBio OTUs
awk '{if($3 == 100) printf $1"\n"}' $OUT > $LIST
sequenceSelect.py -f "pr2_otus99.fasta" -l $LIST -a r -o "pr2_otus99_clean.fasta"

# Add an identifier for future references to both datasets
sed -i "s/^>/>PacBio_/g" "pacbio.fasta"
sed -i "s/^>/>PR2_/g" "pr2_otus99_clean.fasta"

# And finally merge
cat "pacbio.fasta" "pr2_otus99_clean.fasta" > pb_pr2.fasta
