#!/bin/bash

# Set file names -----------------------------------------------------------------------------------
# The directory where the trees and control files are stored
toDATE="toDateD/"
# The prefix of the tree and control files to be dated
FILE="step3f_RAcat1_rep1_rootD_renamed"
# A directory to store the phylogenetic tree dated
OUTDIR="RAcat11"
# The number of dating replicates
ITERATIONS="100"

TREEFILE="$FILE.tre"
CONTROL="$FILE.ctl"
OUTFILE=${TREEFILE/.tre/_treePLdated.tre}

# Copy the corresponding tree and control file in the new directory --------------------------------
[ ! -d $OUTDIR ] && mkdir -p $OUTDIR
cp "$toDATE/$FILE*" "$OUTDIR/."

cd $OUTDIR

# Insert the tree file name in the option of the control file --------------------------------------
sed -i "s/treefile.*/treefile = $TREEFILE/g" $CONTROL

# And run the different replicates for the given tree ----------------------------------------------
for ((i=1;i<=$ITERATIONS;i++)); do
	echo -e "  Replicate $i\n"
	OUTFILEi=${OUTFILE/.tre/$i.tre}
	sed -i "s\outfile.*\outfile = $OUTFILEi\g" $CONTROL
	treePL $CONTROL
	cat $OUTFILEi >> $OUTFILE
	rm -f $OUTFILEi ${OUTFILEi/.tre/.tre.r8s}
	echo ""
done

# Change the last output file name for future references -------------------------------------------
sed -i "s\outfile.*\outfile = $OUTFILE\g" $CONTROL

# Transform the file into nexus to be read by treeAnnotator ----------------------------------------
echo -e "  Transforming into Nexus format\n"
newick2nexus.py -t "$OUTFILE" -o "${OUTFILE/.tre/.nex}"

# Summarise all trees into one tree with the median node heights -----------------------------------
echo -e "  Annotating nodes\n"
treeannotator -heights median "${OUTFILE/.tre/.nex}" "${OUTFILE/.tre/_median.tre}"

# Set file names -----------------------------------------------------------------------------------
echo "Date trees written to: $OUTFILE"
echo "Final tree written to: ${OUTFILE/.tre/_median.tre}"
echo "Summary control file written to: $CONTROL"
