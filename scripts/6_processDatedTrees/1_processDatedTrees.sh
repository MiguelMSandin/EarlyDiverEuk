#!/bin/bash

# The dated and summarised tree --------------------------------------------------------------------
TREE=$(ls *median.tre)
# The reverse list of IDs to reverse
IDsREV="../../tipNamesIDs_reverse.tsv"

# Copy the resulting tree to keep a raw copy -------------------------------------------------------
cp $TREE ${TREE/.tre/_raw.tre}

# --------------------------------------------------------------------------------------------------
# ------------------------ Annotate manually the tree and save a newick file -----------------------
# ----- Please note that this script assumes the newick file is stored as ${TREE/.tre/.newick} -----
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# -------------------- In this repository all trees are in one common directory --------------------
# -- I recomend you to create one directory for each tree since there will be many files per tree --
# --------------------------------------------------------------------------------------------------

# fix the tree if needed
fixTree.sh $TREE

# Run treeLTTsubTrees.R to get LTT data for every subclade -----------------------------------------
TITLE=$(pwd | sed 's\.*/\\g')
treeLTTsubTrees.R -t $TREE -p ${TREE/.tre/_LTT-subtrees.pdf} -c 300 -l eukaryotes -k eukaryotes -T $TITLE

# Check the median time (and max and min HPD) for the root of the tree -----------------------------
head -n 2 ${TREE/.tre/_LTT-data.tsv} | cut -f 2,4,5

# Extract subclades to a directory -----------------------------------------------------------------
treeExtractSubtrees.R -t $TREE -m 300 -d "clades"

# Change name of the tree tips if you want to check for specific tips/clades/nodes -----------------
python3.12 /usr/local/bin/treeTipRename.py -t ${TREE/.tre/.newick} -l $IDsREV

# Colour the newick tree ---------------------------------------------------------------------------
python3.12 /usr/local/bin/treeColourBranches.py -t ${TREE/.tre/_renamed.newick} -c eukprot -k

# And keep the directory cleaned -------------------------------------------------------------------
rm -f ${TREE/.tre/_renamed.newick}
mv ${TREE/.tre/_renamed_coloured.newick} ${TREE/.tre/_figure.tre}
