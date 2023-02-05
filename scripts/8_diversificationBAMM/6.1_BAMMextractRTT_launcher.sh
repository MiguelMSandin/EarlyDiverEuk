#!/bin/bash

# Here we will substitute the dated subclade tree
TREE="INSERT_TREE_HERE"
# Here we will substitute the given *event_data* file
EDATA="INSERT_EDATA_HERE"
# And here we extract the given number of shifts, that should be in the eventdata file
SHIFT=$(echo $EDATA | sed 's/.*shifts//g' | sed 's/_.*//g')

# And we run the Rscript to extract the Rate Through Time plot and best shift configuration dates
Rscript "PATH/TO/SCRIPT/6.1_rttBAMMextract.R" -t $TREE -e $EDATA -s $SHIFT
