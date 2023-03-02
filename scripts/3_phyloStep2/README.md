# 3_phyloStep2/

This step and the following ([3_phyloStep2](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/3_phyloStep2)) are aimed at completing the phylogenetic trees based on biological relevance. Here we will include those OTUs representatives of at least 2 other sequences or with at least 2 reads. And given the relatively large number of sequences we will use the '*-fast*' option of IQ-Tree to infer the phylogenetic trees.  

[0_alignAndTrim.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/3_phyloStep2/0_alignAndTrim.sh):  
First, we will align and trim both datasets, the forward and the reverse.  

[1_phyloIQTree.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/3_phyloStep2/1_phyloIQTree.sh):  
Then, we will infer the phylogenetic trees under the GTR + Gamma + Invariant sites (```-m GTR+I+G```) using as topological constraint those trees obtained in the first step and with 2 replicates.   

[2_processTrees.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/3_phyloStep2/2_processTrees.sh):  
And lastly, we will remove long branches (with the gESD method as described in the first step) and intruders (those sequences classified in a given supergroup that are resolved within another supergroup).  
