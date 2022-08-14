#!/bin/bash

TREE=""
DIR=""

# First, find the logLikelihood of the tree --------------------------------------------------------
grep "SCORE" ${TREE/.tre/.log}
grep "Total CPU" ${TREE/.tre/.log}
grep "required" ${TREE/.tre/.log}

# Then, prune the tree from long branches ----------------------------------------------------------
treePruneOutliers.py -t $TREE -m "gesd" -i -T "10"

# Root the tree so we can check for inruders -------------------------------------------------------
treeRootOutgroup.py -t ${TREE/.tre/_pruned-gesd.tre} -p Discoba

treeCheckIntruders.py -t ${TREE/.tre/_pruned-gesd_rooted.tre} -a supergroup
# Check the intruders, and add tips if needed ------------------------------------------------------
# If more tip names were added, run from here -----------------------------------------------------
treePruneList.py -t ${TREE/.tre/_pruned-gesd.tre} -l ${TREE/.tre/_pruned-gesd_rooted_intruders.list}

# Rename the final tree to ease the reading of the files -------------------------------------------
mv ${TREE/.tre/_pruned-gesd_pruned.tre} ${TREE/.tre/_cleaned.tre}

# Root again the final tree with Discoba and Amorphea ----------------------------------------------
treeRootOutgroup.py -t ${TREE/.tre/_cleaned.tre} -p Discoba
mv ${TREE/.tre/_cleaned_rooted.tre} ${TREE/.tre/_cleaned_rootD.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootD.tre} -c eukprot

treeRootOutgroup.py -t ${TREE/.tre/_cleaned_rootD.tre} -p Opisthokonta Amoebozoa
mv ${TREE/.tre/_cleaned_rootD_rooted.tre} ${TREE/.tre/_cleaned_rootA.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootA.tre} -c eukprot

[ ! -d "$DIR" ] && mkdir -p "$DIR"
mv ${TREE/.tre/_pruned*} $DIR/
