#!/bin/bash

cd "MC01r"

TREEs="../trees/rootAr/"
CONTROL="pb_pr2_step3_treePL_dating.ctl"
toDATE="toDate-A/"
IDs="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNamesIDs2.tsv"

[ ! -d "$toDATE" ] && mkdir -p "$toDATE"

# Rename the trees to safe memory
for FILE in $(ls $TREEs); do
	echo "Renaming $FILE"
	treeTipRename.py -t "$TREEs$FILE" -l $IDs
	mv $TREEs${FILE/.tre/_renamed.tre} "$toDATE/${FILE/.tre/_renamed.tre}"
	echo ""
done

# Simplify the names of the files
for FILE in $(ls $toDATE); do mv "$toDATE$FILE" "$toDATE${FILE/iqtreef_GTRg_/}"; done
for FILE in $(ls $toDATE); do mv "$toDATE$FILE" "$toDATE${FILE/_cleaned/}"; done

# And now rename the control file to match the renamed tree tips
cp ../$CONTROL .
fileNameReplace.py -f $CONTROL -l $IDs

# Check if all renaming worked
treePL_testControlFile.py -c ${CONTROL/.ctl/_renamed.ctl} -t "$toDATE"$(ls "$toDATE" | head -n 1) -e

# Create a control file for every tree
for FILE in $(ls $toDATE*renamed.tre); do
	echo $FILE
	TREE=$(basename -- $FILE)
	OUT=${TREE/.tre/_dated.tre}
	cp ${CONTROL/.ctl/_renamed.ctl} $toDATE${TREE/.tre/.ctl}
	sed -i "s/treefile =/treefile = $TREE/g" $toDATE${TREE/.tre/.ctl}
	sed -i "s/outfile =/outfile = $OUT/g" $toDATE${TREE/.tre/.ctl}
done

# And start dating...
# cd toDate
# for FILE in $(find *ctl); do
# 	echo "Dating ${FILE/.ctl/.tre}"
# 	treePL $FILE
# 	rm -f *r8s
# 	mv ${FILE/.ctl/_dated.tre} ../dated/${FILE/.ctl/_dated.tre}
# done

# Or better use the "treeStepDating_3_treePL_replicates.sh" to run 100 replicates of each tree
