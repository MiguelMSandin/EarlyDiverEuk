#!/bin/bash

cd "MC01r/rootD/XXXX"

TREE=""

cp $TREE ${TREE/.tre/_raw.tre}

IDs="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNamesIDs.tsv"
IDsREV="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNamesIDs_reverse.tsv"

# Annotate manually the tree, save a copy and save a newick file
fixTreeCommented.sh $TREE

# Run LTT-subTrees.R to get LTT data for every subclade
TITLE=$(pwd | sed 's\.*/\\g')
treeLTTsubTrees.R -t $TREE -p ${TREE/.tre/_LTT-subtrees.pdf} -c 300 -l eukaryotes -k eukaryotes -T $TITLE

head -n 2 ${TREE/.tre/_LTT-data.tsv} | cut -f 2,4,5

# Change name of the tree tips
treeTipRename.py -t ${TREE/.tre/.newick} -l $IDsREV

# Colour the newick tree
treeColourBranches.py -t ${TREE/.tre/_renamed.newick} -c eukprot -k

# And keep the directory cleaned
rm -f ${TREE/.tre/_renamed.newick}
mv ${TREE/.tre/_renamed_coloured.newick} ${TREE/.tre/_figure.tre}

# Once everything is finished, extract the root ages for further analysis
# cd ..
# for FILE in $(find */*data.tsv); do
# 	echo "$FILE" $(head -n 2 $FILE | cut -f 2,4,5 | tail -n 1) >> tmp.tsv
# done
