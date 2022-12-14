#!/bin/bash

cd /home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/

TREES="trees/rootDr/"
COMMON="common_r/common_tips.list"
CONTROL="pb_pr2_step3_treePL_dating.ctl"

# Get common tips in all trees
treeCompareTips.py -t $TREES*tre -e "common" -o $COMMON

# Prune all trees so the all have common tips
for TMP in $(ls $TREES); do
	echo $TREES$TMP
	treePruneList.py -t $TREES$TMP -l $COMMON -i -o "common_r/trees/${TMP/.tre/_pruned.tre}"
	treeColourBranches.py -t "common_r/trees/${TMP/.tre/_pruned.tre}" -c eukprot -k
	mv "common_r/trees/${TMP/.tre/_pruned_coloured.tre}" "common_r/trees/coloured/."
	echo ""
done

# Extract number of reads for every sequence
tmp=$(mktemp --tmpdir=$(pwd))
LEN=$(wc -l < $COMMON)
i=0
while read LINE; do
	((i=i+1))
	echo -n -e "\r    $i/$LEN"
	if echo $LINE | grep -q "^PacBio"; then
		echo $LINE | sed 's/.*_Otu[0-9]*_//g' | sed 's/_.*//g' >> $tmp
	elif echo $LINE | grep -q "^PR2";then
		echo $LINE | sed 's/.*_//g' >> $tmp
	fi
done < "${COMMON}"
echo ""
paste $COMMON $tmp > ${COMMON/.list/.tsv}
rm -f $tmp

# Subset only sequences with 10 or more reads to select the most reliable sequences
awk '$2>=10' ${COMMON/.list/.tsv} > ${COMMON/.list/_reads10.tsv}

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ------------------------------- Build the control file carefully ! -------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Quickly check selected tips in all trees
treePL_summaryControlFile.py -t common_r/trees/*tre -c $CONTROL -o ${CONTROL/.ctl/_reverse_summary.tsv}

NAMES=""
treeLCAcount.py -t common_r/trees/*tre -n $NAMES -l tmp.tsv

# Check if the control file is compatible with the trees
for TMP in $(ls common_r/trees/*tre); do
	echo $TMP
	treePL_testControlFile.py -c $CONTROL -t $TMP
	mv ${TMP/.tre/_calibrated.tre} common_r/trees_calibrated/
	echo ""
done
