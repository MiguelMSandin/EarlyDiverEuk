#!/bin/bash

# Here we will substitute the dated subclade tree
TREE="INSERT_TREE_HERE"
# Here we will substitute the given fraction file
FRACTION="INSERT_FRAC_HERE"
# And here we will round up the fraction number for file naming only
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')

# And finally we run the script '1.2_clads.jl'
julia 1.2_clads.jl -t $TREE -f $FRACTION -r ${TREE/.tre/_div${FRAC}_clads.Rdata}
