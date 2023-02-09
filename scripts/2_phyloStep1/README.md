# 2_phyloStep1/

In this step we are going to build the backbone phylogenetic tree including only OTUs representative of at least 10 other sequences. In order to allow phylogenetic uncertainty, we will perform several phylogenetic analyses with different softwares (RAxML and RAxML-ng) and different models of evolution (GTR+Gamma and GTR+CAT).  

[0_alignAndTrim.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/2_phyloStep1/0_alignAndTrim.sh):  
In addition, we will align the sequences over the forward and reverse dataset in parallel to grant the possibility of misalignments of highly divergent sequences.  

[1_phyloRAxML-GTR-GAMMA.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/1_phyloRAxML-GTR-GAMMA.sh):  
We will first perform a thorough search for the maximum likelihood tree over 100 BootStrap replicates in ```RAxML``` including the option ```-D```, which optimizes the ML search convergence criteria for very large datasets and the model of nucleotide substitution GTR+Gamma (```-m GTRGAMMA```).  

[2_phyloRAxML-GTR-CAT.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/2_phyloRAxML-GTR-CAT.sh):  
We will also perform a second thorough search (as above) with the model of nucleotide substitution GTR+CAT ()```-m GTRCAT```).  

[3_phyloRAxML-ng.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/3_phyloRAxML-ng.sh):  
Lastly, we will also perform 10 quick searches in ```RAxML-ng``` and select the 4 best scoring trees to allow a greater phylogenetic uncertainty in deep nodes.  

We will perform the same 6 phylognetic analyses over the reverse alignment:  
[4_phyloRAxML-GTR-GAMMA_reverse.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/4_phyloRAxML-GTR-GAMMA_reverse.sh).  
[5_phyloRAxML-GTR-CAT_reverse.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/5_phyloRAxML-GTR-CAT_reverse.sh).  
[6_phyloRAxML-ng_reverse.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/6_phyloRAxML-ng_reverse.sh).  

[7_processTrees.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/2_phyloStep1/7_processTrees.sh):  
And finally, we will remove long-branches by identifying outliers from a normal distribution (applying the generalized Extreme Studentized Deviate, gESD, method from [Rosner, 1983](https://www.tandfonline.com/doi/abs/10.1080/00401706.1983.10487848)). In addition, we will root the trees and colour the branches for visualization and identification of errors purposes at this state.  
