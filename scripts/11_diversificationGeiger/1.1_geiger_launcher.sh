#!/bin/bash

TREE="INSERT_TREE"
FRACTION=INSERT_FRACTION
CLADE=INSERT_CLADE
THREADS=4
SCRIPT="/home/miguel/Documents/Uppsala/1_ecoEvo/repository/scripts/11_diversificationGeiger/1.2_geiger.R"

FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
PREFIX="${CLADE}_d${FRAC}"

Rscript $SCRIPT -t $TREE -d $FRACTION -n 0.009 -e 0.3 -c $THREADS -p $PREFIX

gzip ${PREFIX}*
