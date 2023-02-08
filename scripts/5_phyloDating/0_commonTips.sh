#!/bin/bash

# We build the control file based on a tree containing tips common to all trees --------------------
# In other words, we remove those tips that do not appear in all trees -----------------------------

# The last directory where you have all phylogenetic trees with the same root
TREES="trees/rootD/"
# A directory to store the common tip names and trees
COMMON="common/"
[ ! -d "$COMMON" ] && mkdir -p "$COMMON"
[ ! -d "$COMMON/trees" ] && mkdir -p "$COMMON/trees"

# Get common tips in all trees
treeCompareTips.py -t $TREES*tre -e "common" -o "$COMMON/common_names.list"

# Prune all trees so they all have common tips
for TMP in $(ls $TREES); do
	echo $TREES$TMP
	treePruneList.py -t $TREES$TMP -l "$COMMON/common_names.list" -i -o "$COMMON/trees/${TMP/.tre/_pruned.tre}"
	treeColourBranches.py -t "$COMMON/trees/${TMP/.tre/_pruned.tre}" -c eukprot -k
	mv "$COMMON/trees/${TMP/.tre/_pruned_coloured.tre}" "$COMMON/trees/coloured/."
	echo ""
done

# Extract number of reads for every sequence
tmp=$(mktemp --tmpdir=$(pwd))
LEN=$(wc -l < "$COMMON/common_names.list")
i=0
while read LINE; do
	((i=i+1))
	echo -n -e "\r    $i/$LEN"
	if echo $LINE | grep -q "^PacBio"; then
		echo $LINE | sed 's/.*_Otu[0-9]*_//g' | sed 's/_.*//g' >> $tmp
	elif echo $LINE | grep -q "^PR2";then
		echo $LINE | sed 's/.*_//g' >> $tmp
	fi
done < "$COMMON/common_names.list"
echo ""
paste "$COMMON/common_names.list" $tmp > "$COMMON/common_names_reads.tsv"
rm -f $tmp

# Subset only sequences with 10 or more reads to select the most reliable sequences for the dating
awk '$2>=10' "$COMMON/common_names_reads.tsv" > "$COMMON/common_names_reads10.tsv}"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ------------------------------- Build the control file carefully ! -------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
