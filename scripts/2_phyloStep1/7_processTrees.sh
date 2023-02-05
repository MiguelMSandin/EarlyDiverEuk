#!/bin/bash

TREE=""
DIR=""

# Prune the tree from long branches ----------------------------------------------------------------
echo -e "\nPruning outliers"
treePruneOutliers.py -t $TREE -m "gesd" -i -T "10"

# There should be no intruders, so there is no need to check in advanced
# Prune the unrooted and pruned tree for downstream analysis ---------------------------------------
echo -e "\nPruning intruders"
treePruneList.py -t ${TREE/.tre/_pruned-gesd.tre} -l ${TREE/.tre/_pruned-gesd_rooted_intruders.list}

# Rename the final tree to ease the reading of the files -------------------------------------------
echo -e "\nRenaming '${TREE/.tre/_pruned-gesd_pruned.tre}' to '${TREE/.tre/_cleaned.tre}'"
mv ${TREE/.tre/_pruned-gesd_pruned.tre} ${TREE/.tre/_cleaned.tre}

# Root again the final tree in Discoba -------------------------------------------------------------
echo -e "\nRooting in Discoba for visualization purposes"
treeRootOutgroup.py -t ${TREE/.tre/_cleaned.tre} -p Discoba
mv ${TREE/.tre/_cleaned_rooted.tre} ${TREE/.tre/_cleaned_rootD.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootD.tre} -c eukprot

# And now root in Amorphea
echo -e "\nRooting in Amorphea for visualization purposes"
treeRootOutgroup.py -t ${TREE/.tre/_cleaned_rootD.tre} -p Opisthokonta Amoebozoa
mv ${TREE/.tre/_cleaned_rootD_rooted.tre} ${TREE/.tre/_cleaned_rootA.tre}
treeColourBranches.py -t ${TREE/.tre/_cleaned_rootA.tre} -c eukprot

# Check the final tree -----------------------------------------------------------------------------
# And moving intermediate files for clarity --------------------------------------------------------
mkdir $DIR
mv *$DIR*pruned* $DIR/
treeRemoveBranchLengths.py -t ${TREE/.tre/_cleaned.tre}
