#!/bin/bash

# The dated and summarised clade -------------------------------------------------------------------
CLADE="clade_Archaeplastida.tre"

cp $CLADE ${CLADE/.tre/_annot.tre}
# python3.12 /usr/local/bin/treeAnnotate.py -t clade_Archaeplastida.tre -n Chlorophyta Streptophyta Rhodophyta Picozoa Glaucophyta -e 0.95

# --------------------------------------------------------------------------------------------------
# ------------------------ Annotate manually the tree and save a newick file -----------------------
# ----- Please note that this script assumes the newick file is stored as ${TREE/.tre/.newick} -----
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# -------------------- In this repository all trees are in one common directory --------------------
# -- I recomend you to create one directory for each tree since there will be many files per tree --
# --------------------------------------------------------------------------------------------------

# fix the tree if needed
fixTree.sh ${CLADE/.tre/_annot.tre}

# Run treeLTTsubTrees.R to get LTT data for every subclade -----------------------------------------
treeLTT.py -t ${CLADE/.tre/_annot.tre} -f nexus -s
