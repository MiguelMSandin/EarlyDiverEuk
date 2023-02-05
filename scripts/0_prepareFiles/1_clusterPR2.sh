#! bin/bash

THREADS="8"

# First extract sequences for every supergroup in order to decrease RAM memory usage
# Still, this might be relatively demanding in terms of RAM memory for certain groups
# I would suggest split the loop in two (first extract groups, then align), and ask for appropriate RAM memory according to the size of the fasta files
for GROUP in $(cut -f 3 -d '_' pr2_raw_names.tsv | sort | uniq); do
	FASTA="subset_"$GROUP".fasta"
	ALIGN="subset_"$GROUP"_align.fasta"
	OUT="otus99_"$GROUP".fasta"
	sequenceSelect.py -f "pr2_raw.fasta" -p $GROUP -o "subset_$GROUP.fasta" -v
	# Align the sequences
	mafft --thread $THREADS $FASTA > $ALIGN
	# Calculate distance
	mothur "#dist.seqs(fasta=$ALIGN, cutoff=0.1, output=lt, processors=$THREADS)"
	# Cluster sequences at 99% similariy (or 1% dissimilarity)
	mothur "#cluster(phylip=${ALIGN/.fasta/.phylip.dist}, cutoff=0.01)"
	# And get the centroid of each cluster as representative sequence
	mothur "#get.oturep(phylip=${ALIGN/.fasta/.phylip.dist}, list=${ALIGN/.fasta/.phylip.opti_mcc.list}, fasta=$FASTA, cutoff=0.01, method=distance)"
	# Change name of final OTU file of the group
	mv ${ALIGN/.fasta/.phylip.opti_mcc.0.01.rep.fasta} $OUT
	# And delete temporary files
	rm -f mothur* subset*
done

# cat all otu files into one single file
cat "otus99*" > "pr2_otus99.fasta"

# And replace the tabs separating the OTU name and the pipe
sed -i "s/\t/_/g" "pr2_otus99.fasta"
tmp1=$(mktemp --tmpdir=$(pwd))
cat "pr2_otus99.fasta" | tr '|' '_' > $tmp1
mv $tmp1 "pr2_otus99.fasta"
