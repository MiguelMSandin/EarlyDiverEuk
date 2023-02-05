#!/bin/bash

# cd to the root directory of a dated phylogenetic tree, where you have, at least, the subclades and the diversity fractions file
# First create a folder for CLaDS analyses
DIR="clads"
[ ! -d "$DIR" ] && mkdir -p "$DIR"
# A prefix to name two tables for the minimum ("$PREFIX"_min.tsv) and maximum ("$PREFIX"_max.tsv) diversity fractions for further analyses
PREFIX="fractions"

# The diversity fractions file
FRACTIONS_RAW=$(find fractions/*fractions.tsv)
# All dated trees of the subclades
TREES=$(ls clades/clade*)
# And get all clade names from the trees
CLADES=$(echo $TREES | sed 's\clades/clade_\\g' | sed 's/\.tre//g')

# Loop through all subclades to extract minimum and maximum diversity fractions
for CLADE in $CLADES; do
	min=$(grep $CLADE $FRACTIONS_RAW | cut -f 4 | sort -n | head -n 1)
	max=$(grep $CLADE $FRACTIONS_RAW | cut -f 4 | sort -n | tail -n 1)
	TREE="../clades/"$(ls clades | grep $CLADE)
	echo -e "$CLADE\t$TREE\t$min" >> $DIR/"$PREFIX"_min.tsv
	echo -e "$CLADE\t$TREE\t$max" >> $DIR/"$PREFIX"_max.tsv
done

# Finally copy all subclade trees to the '$DIR' directory to prepare for further analyses
cp clades/clade_* clads/.
