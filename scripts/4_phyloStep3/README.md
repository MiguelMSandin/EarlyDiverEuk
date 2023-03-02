# 4_phyloStep3/

This step is aimed at completing the phylogenetic trees created until now with all the rest of the OTUS, which are those representatives of 1 sequences or with at least 1 read. And given the relatively large number of sequences we will use the ```-fast``` option of IQ-Tree to infer the phylogenetic trees, as in [3_phyloStep2](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/3_phyloStep2).  

As in [3_phyloStep2](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/3_phyloStep2), we will align, trim, infer the phylogenetic trees and process them using the same protocol:  
[0_alignAndTrim.sh](0_alignAndTrim.sh):  
[1_phyloIQTree.sh](1_phyloIQTree.sh):  
[2_processTrees.sh](2_processTrees.sh):  
  
[3_infoTrees.R](3_infoTrees.R):  
After completing all phylogenetic analyses, we will quickly explore the likelihood of all trees and the variance in the number of OTUs per supergroup.  
