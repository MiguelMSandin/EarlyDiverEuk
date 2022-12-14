#!/bin/bash

TREE=""
DIR=""

# First, find the logLikelihood of the tree --------------------------------------------------------
echo -e "\nLikelihood: $(grep "SCORE" ${TREE/.tre/.log})"
echo "$(grep "Total CPU" ${TREE/.tre/.log})"
echo "RAM memory required: $(grep "required" ${TREE/.tre/.log})"

# Then, prune the tree from long branches ----------------------------------------------------------
echo -e "\nPruning outliers"
treePruneOutliers.py -t $TREE -m "gesd" -i -T "10"

# Root the tree so we can check for inruders -------------------------------------------------------
echo -e "\nChecking for intruders (after rooting in Discoba)"
treeRootOutgroup.py -t ${TREE/.tre/_pruned-gesd.tre} -p Discoba
treeCheckIntruders.py -t ${TREE/.tre/_pruned-gesd_rooted.tre} -a supergroup

# Prune the unrooted and pruned tree for downstream analysis ---------------------------------------
# If you have added tip names to the intruderes name list, run from here ---------------------------
echo -e "\nPruning intruders"
treePruneList.py -t ${TREE/.tre/_pruned-gesd.tre} -l ${TREE/.tre/_pruned-gesd_rooted_intruders.list}

# Rename the final tree to ease the reading of the files -------------------------------------------
echo -e "\nRenaming '${TREE/.tre/_pruned-gesd_pruned.tre}' to '${TREE/.tre/_cleaned.tre}'"
mv ${TREE/.tre/_pruned-gesd_pruned.tre} ${TREE/.tre/_cleaned.tre}

# Root again the final tree with Discoba and Amorphea ----------------------------------------------
echo -e "\nRooting in Discoba for visualization purposes"
treeRootOutgroup.py -t ${TREE/.tre/_cleaned.tre} -p Discoba
mv ${TREE/.tre/_cleaned_rooted.tre} ${TREE/.tre/_cleaned_rootD.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootD.tre} -c eukprot

echo -e "\nRooting in Amorphea for visualization purposes"
treeRootOutgroup.py -t ${TREE/.tre/_cleaned_rootD.tre} -p Opisthokonta Amoebozoa
mv ${TREE/.tre/_cleaned_rootD_rooted.tre} ${TREE/.tre/_cleaned_rootA.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootA.tre} -c eukprot

# And moving intermediate files for clarity --------------------------------------------------------
mkdir $DIR
mv *$DIR*pruned* $DIR/
treeRemoveBranchLengths.py -t ${TREE/.tre/_cleaned.tre}

echo -e "\nIntruders: $(wc -l $DIR/${TREE/.tre/_pruned-gesd_rooted_intruders.list})\n"
echo -e "Tips:\n$(treeCountTips.py -t $TREE ${TREE/.tre/_cleaned.tre})"


# fixTree.sh ${TREE/.tre/_cleaned_rootD_coloured.tre}
# fixTree.sh ${TREE/.tre/_cleaned_rootA_coloured.tre}
