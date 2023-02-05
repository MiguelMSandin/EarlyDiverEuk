#!/bin/bash

# 'cd' to the 'clads' directory created with the previous script

# This is the script we will call to extract the data from the Rdata object from the CLaDS analyses
SCRIPT="PATH/TO/SCRIPT/2.1_rttClaDSextract.R"

# Now will loop through all the *Rdata files to extract the RTT and DTT data, create a directory and move the file
for FILE in $(ls | grep Rdata); do
	CLADE=$(echo $FILE | sed 's/clade_//g' | sed 's/_.*//g')
	# If you want to run more analysis, comment the next line so the subclade tree is not removed
	rm -f "clade_$CLADE.tre"
	[ ! -d "$CLADE" ] && mkdir -p "$CLADE"
	Rscript "$SCRIPT" -f $FILE
	mv ${FILE/.Rdata/*} $CLADE
done
