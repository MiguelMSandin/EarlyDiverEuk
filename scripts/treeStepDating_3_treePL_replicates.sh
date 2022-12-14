#!/bin/bash

FILE="step3f_RAcat1_rep1_rootD_renamed"
OUTDIR="RAcat11"
ITERATIONS="100"

TREEFILE="$FILE.tre"
CONTROL="$FILE.ctl"
OUTFILE=${TREEFILE/.tre/_treePLdated.tre}

[ ! -d $OUTDIR ] && mkdir -p $OUTDIR

echo -e "  Copying files to folder\n"
cp $CONTROL $OUTDIR/$CONTROL
cp $TREEFILE $OUTDIR/$TREEFILE

cd $OUTDIR

sed -i "s/treefile.*/treefile = $TREEFILE/g" $CONTROL

echo -e "  Running treePL with $ITERATIONS replicates\n"

for ((i=1;i<=$ITERATIONS;i++));
do
	echo -e "  Replicate $i\n"
	
	OUTFILEi=${OUTFILE/.tre/$i.tre}
	sed -i "s\outfile.*\outfile = $OUTFILEi\g" $CONTROL
	
	treePL $CONTROL
	
	cat $OUTFILEi >> $OUTFILE
	
	rm -f $OUTFILEi ${OUTFILEi/.tre/.tre.r8s}
	
	echo ""
done

sed -i "s\outfile.*\outfile = $OUTFILE\g" $CONTROL

echo -e "  Transforming into Nexus format\n"
newick2nexus.py -t "$OUTFILE" -o "${OUTFILE/.tre/.nex}"

echo -e "  Annotating nodes\n"
/home/miguel/softwares/BEAST/BEASTv1.10.4/bin/treeannotator -heights median "${OUTFILE/.tre/.nex}" "${OUTFILE/.tre/_median.tre}"
echo -e "\n"
echo "Date trees written to: $OUTFILE"
echo "Final tree written to: ${OUTFILE/.tre/_median.tre}"
echo "Summary control file written to: $CONTROL"
echo -e "\nDone"
