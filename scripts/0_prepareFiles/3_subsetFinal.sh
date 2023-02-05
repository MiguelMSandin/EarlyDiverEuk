#!/bin/bash

# First, extract number of reads of PR2 sequences
tmp1=$(mktemp --tmpdir=$(pwd))
grep "^>PR2_" "pb_pr2.fasta" | sed "s/^>//g" > $tmp1
tmp2=$(mktemp --tmpdir=$(pwd))
sed "s/.*_Otu[0-9]\+_//g" $tmp1 > $tmp2

tmp3=$(mktemp --tmpdir=$(pwd))
paste $tmp1 $tmp2 > $tmp3

# Now of PacBio sequences
grep "^>PacBio_" "pb_pr2.fasta" | sed "s/^>//g" > $tmp1
sed "s/.*_Otu[0-9]\+_//g" $tmp1 | sed "s/_.*//g" > $tmp2

tmp4=$(mktemp --tmpdir=$(pwd))
paste $tmp1 $tmp2 > $tmp4

cat $tmp3 $tmp4 > "otu_reads.tsv"
rm -f tmp*

# Now get those OTUs with 10 or more reads, complemented with clades that do not have OTU such abundant
sequenceSelect.py -f "pb_pr2.fasta" -l "otus_reads10_complemented.list" -o "pb_pr2_reads10c.fasta"

# Now get those OTUs with 10 or more reads
sequenceSelect.py -f "pb_pr2.fasta" -l "otus_reads2_complemented.list" -o "pb_pr2_reads2c.fasta"

# And finally organize the data into folders
mkdir ../
