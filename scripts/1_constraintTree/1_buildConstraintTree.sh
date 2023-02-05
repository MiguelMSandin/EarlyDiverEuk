#!/bin/bash

# Download in-house scripts ------------------------------------------------------------------------
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/checkConstrainTree.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/checkConstrainTaxa.sh
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/buildConstrainTree.sh
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/findSeqs.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/treeStats.py

chmod +x *py
mv *py ../

# Set file names -----------------------------------------------------------------------------------
FASTA="pb_pr2_reads10c.fasta"
TREE="eukProt_v9.tre"
CONSTRAINT_TREE="eukProt_v9_forConstraint.tre"

# Check if the reference tree can be used for the given fasta file _________________________________
checkConstrainTree.py -f $FASTA -t $TREE -d "taxa"
# Check manually every file for correct sequences retrieved
# grep -v Oligohymenophorea taxa/Tardigrada.txt > tmp; mv tmp taxa/Tardigrada.txt
# grep -v 'MAST-1[0-9]' taxa/MAST-1.txt > tmp; mv tmp taxa/MAST-1.txt

# Check if all the taxa we have previously retrieved, is monophyletic and there are no incongruencies
# Be aware that this script is highly dependent of your sequence names. It might not be useful at all!!
checkConstrainTaxa.sh -d "taxa" -o ${CONSTRAINT_TREE/.tre/_checking.log}

# Build the constraint tree ________________________________________________________________________
buildConstrainTree.sh -d "taxa" -t $TREE -o $CONSTRAINT_TREE -n
buildConstrainTree.sh -d "taxa" -t $TREE -o ${CONSTRAINT_TREE/.tre/_annot.tre}

# Check that you have all tips and sequences in the other file _____________________________________
findSeqs.py -f $FASTA -t $CONSTRAINT_TREE -c t
# If all tips are found in the fasta file, then GOOD TO GO!

# Check how many tips and nodes there are
treeStats.py -t $CONSTRAINT_TREE
