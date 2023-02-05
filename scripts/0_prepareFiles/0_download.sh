#! bin/bash

# Download in-house scripts ------------------------------------------------------------------------
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/fasta/fastaRename.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/fasta/sequenceSelect.py

chmod +x *py
mv *py ../

# Download PR2 database ----------------------------------------------------------------------------
wget https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_taxo_long.fasta.gz
gunzip -k "pr2_version_4.14.0_SSU_taxo_long.fasta.gz"
mv "pr2_version_4.14.0_SSU_taxo_long.fasta" "pr2_raw.fasta"
sed -i 's/ /_/g' "pr2_raw.fasta"

# Extract only nuclear genes
tmp1=$(mktemp --tmpdir=$(pwd))
sequenceSelect.py -f "pr2_raw.fasta" -p nucleus -o $tmp1
cp $tmp1 "pr2_raw.fasta"

# Modify names
# First extract the sequence identifiers
grep "^>" "pr2_raw.fasta" | sed 's/>//g' > "pr2_names.tsv"
# Remove extra information from the accession number
cut -f 1 -d '|' "pr2_names.tsv" | sed 's/\..*//g' > "$tmp1"
# Get only the taxonomy columns
tmp2=$(mktemp --tmpdir=$(pwd))
cut -f 6-12 -d '|' "pr2_names.tsv" | tr '|' '_' > "$tmp2"
# Paste the new identifiers
tmp3=$(mktemp --tmpdir=$(pwd))
paste "$tmp1" "$tmp2" -d '_' > "$tmp3"
# Now create a table with the old and new identifiers
paste "pr2_names.tsv" $tmp3 > "$tmp1"
mv "$tmp1" "pr2_names.tsv"
# And change the names
fastaRename.py -f "pr2_raw.fasta" -l "pr2_names.tsv" -o $tmp2
mv "$tmp2" "pr2_raw.fasta"
# Delete dots from the identifiers
sed -i 's/\.//g' "pr2_raw.fasta"
# Remove temporary files
rm -f tmp*
# And finally create a table with the names for future references
grep "^>" "pr2_raw.fasta" | sed 's/>//g' > "pr2_raw_names.tsv"

# Download PacBio database -------------------------------------------------------------------------
# First the 18S
wget https://figshare.com/ndownloader/files/29133525
mv "29133525" "pacbio_18S.fasta.gz"
gunzip "pacbio_18S.fasta.gz"
# Then the 28S
wget https://figshare.com/ndownloader/files/29133528
mv "29133528" "pacbio_28S.fasta.gz"
gunzip "pacbio_28S.fasta.gz"
# Concatenate
fastaConcat.py -f "pacbio_18S.fasta" "pacbio_28S.fasta" -o "pacbio.fasta"

rm -f "pacbio_18S.fasta" "pacbio_28S.fasta"
