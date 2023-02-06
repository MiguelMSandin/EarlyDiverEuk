#!/bin/bash

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------- After building the control file carefully! ---------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# The last directory where you have all phylogenetic trees with the same root
TREES="trees/rootD/"
# The directory where the trees with common tip names are stored
COMMON="common/trees"
# The control file
CONTROL="pb_pr2_treePL_dating.ctl"
[ ! -d "$COMMON/trees_calibrated" ] && mkdir -p "$COMMON/trees_calibrated"

# Quickly check selected tips in all trees
1.1_treePL_summaryControlFile.py -t "$COMMON/*tre" -c $CONTROL -o "${CONTROL/.ctl/_summary.tsv}"
# If you find problems, you can quickly spot them by selecting the tip names in the following script
# NAMES=""
# treeLCAcount.py -t common_r/trees/*tre -n $NAMES -l tmp.tsv

# Check if the control file is compatible with all the trees
for TMP in $(ls "$COMMON/*tre"); do
	echo $TMP
	1.2_treePL_testControlFile.py -c $CONTROL -t $TMP
	mv ${TMP/.tre/_calibrated.tre} "$COMMON/trees_calibrated"
	echo ""
done

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# -------------------------- If everything seems to be fine, ready to go! --------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
