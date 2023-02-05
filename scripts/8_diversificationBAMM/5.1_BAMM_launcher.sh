#!/bin/bash

# Here we will substitute the given control file
FILE="INSERT_FILE_NAME_HERE"
OUT=$(grep outName $FILE | cut -f 3 -d ' ')

# Run BAMM
bamm -c $FILE

# Compress all files to safe memory
gzip $OUT*
# And move them file to a new folder for clarity
mv $OUT* ../../diver/
