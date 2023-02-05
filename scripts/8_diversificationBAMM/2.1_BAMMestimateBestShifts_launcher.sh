#!/bin/bash

# Here we will substitute the given control file
FILE="INSERT_FILE_NAME_HERE"
OUT=$(grep outName $FILE | cut -f 3 -d ' ')

# Run BAMM
bamm -c $FILE

# We remove all files except the mcmc chain and compress the mcmc chain file, to safe memory
rm -f $OUT*chain_swap.txt
rm -f $OUT*event_data.txt
rm -f $OUT*run_info.txt
gzip $OUT*

# And move the output file to a new folder for clarity
mv $OUT* "../../diver/"
