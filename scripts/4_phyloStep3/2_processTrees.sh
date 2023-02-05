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
# treeRootOutgroup.py -t ${TREE/.tre/_pruned-gesd.tre} -l "root_discoba.list"

# Check for intruders ------------------------------------------------------------------------------
treeCheckIntruders.py -t ${TREE/.tre/_pruned-gesd_rooted.tre} -a supergroup
# Check the intruders, and add tips if needed ------------------------------------------------------
# If more tip names were added, run from here ------------------------------------------------------
treePruneList.py -t ${TREE/.tre/_pruned-gesd.tre} -l ${TREE/.tre/_pruned-gesd_rooted_intruders.list}

# Rename the final tree to ease the reading of the files -------------------------------------------
mv ${TREE/.tre/_pruned-gesd_pruned.tre} ${TREE/.tre/_cleaned.tre}

# Root again the final tree with Discoba and colour the branches -----------------------------------
treeRootOutgroup.py -t ${TREE/.tre/_cleaned.tre} -p Discoba
mv ${TREE/.tre/_cleaned_rooted.tre} ${TREE/.tre/_cleaned_rootD.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootD.tre} -c eukprot -k

# And now the same with Amorphea as a root
treeRootOutgroup.py -t ${TREE/.tre/_cleaned_rootD.tre} -p Fungi Metazoa Amoebozoa
mv ${TREE/.tre/_cleaned_rootD_rooted.tre} ${TREE/.tre/_cleaned_rootA.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootA.tre} -c eukprot -k

# And keep the directory clean
[ ! -d "$DIR" ] && mkdir -p "$DIR"
mv ${TREE/.tre/_pruned*} $DIR/
[ ! -d "rootD" ] && mkdir -p "rootD"
mv ${TREE/.tre/_cleaned_rootD*} "rootD/"
[ ! -d "rootA" ] && mkdir -p "rootA"
mv ${TREE/.tre/_cleaned_rootA*} "rootA/"

# And copy the raw trees to the dating folder
cp rootD/${TREE/.tre/_cleaned_rootD.tre} ../../stepDating/trees/rootD
cp rootA/${TREE/.tre/_cleaned_rootA.tre} ../../stepDating/trees/rootA
