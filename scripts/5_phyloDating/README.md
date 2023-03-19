# 5_phyloDating

After resolving all phylogenetic trees, we will begin the dating. Here we will relate the DNA sequence extracted from a living organism with a (most likely) extinct fossil group. In some cases it is possible to find a continuous fossil record with varying but related morphologies, and thus the dating is somehow straight forward. In other cases, such link is not as obvious, and an expertise is needed to determine what extinct morphology corresponds (or rather it is assumed to correspond) to the first representative in the fossil record for a given group. The older we go back in time the more difficult will be to establish such a link.  
In this context, it is very important to do an extensive literature research on what nodes can be calibrated, have been previously used by experts in the group to calibrate equivalent nodes, and make sense to be calibrated from the rDNA phylogenetic point of view.  
Briefly, it is worth dedicating a big effort at this step to make sure the calibration is consistent with the fossil record, the phylogenetic trees and the living diversity such calibration gather. So probably you will need to go over several runs of trial-and-error until you have a consistent calibration.  
Once we have a table with all possible calibrations (in this case: the node and minimum and/or maximum dates of appearance) we can begin with the calibration.  
Given the large number of tips, applying a bayesian approach will be a computational effort that most likely will never arrive to a convergence. Therefore, we will apply a Penalized Likelihood approach implemented in [TreePL](https://github.com/blackrim/treePL), that directly calibrates the phylogenetic distance of the trees. This software requires a control file with the parameters of the run and the calibrations.  
  
[0_commonTips.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/0_commonTips.sh):  
In order to build a solid control file that can be used for all trees (and therefore accounting for phylogenetic uncertainty) we will first prune all tips that do **not** appear in all trees. Using several random pruned trees we will build **carefully** the control file paying special attention to  
- the node is monophyletic,  
- the node is not ambiguous (i.e.; consistent among phylogenetic trees),  
- the calibration node covers at least 4 OTU sequences in all trees (an arbitrary threshold),  
- the branch-lengths are consistent with the calibration (i.e.; if a node gathers 4 sequences with near-0 branch length, a calibration of >100 million years is doubtful; if a given node gathers 20 sequences with very long branches, a calibration of <5 million years is doubtful), **and**  
- the calibration has been previously used and validated in clade-specific analyses.  
  
[1_CheckControl.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/1_CheckControl.sh):  
Here we will check that the control file is OK. We are interested in checking:  
- [1.1_treePL_summaryControlFile.py](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/): that every calibrated node gathers more or less the same number of tips, so it is easier to spot paraphyly or not consistent tips among different phylogenetic trees, and  
- [1.2_treePL_testControlFile.py](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/1.2_treePL_testControlFile.py): that the control file is compatible with all trees. That is that every taxon used to get the last common ancestor of the node is present in all trees and that the calibrations are compatible with one another (i.e.; parent nodes with older -or equal- ages than child nodes).  
  
If everything is OK until this point, we are finally ready to date the trees.  
[2_multiplyControl.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/2_multiplyControl.sh): But first we have to create different control files for each phylogenetic tree. Besides, since the tip names are very long, we will replace them for an identifier to safe memory.  
  
[3_treePL_replicates.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/5_phyloDating/3_treePL_replicates.sh): Finally, we can date the trees. TreePL is a rather fast algorithm, so we will benefit from that and replicate every dating analysis 100 times to allow uncertainty. Then, we will summarise all 100 replicates with [treeannotator](https://beast.community/treeannotator) and use the median as a final node height.  
  
At the end we will have 96 total dated phylogenetic trees of eukaryotes relatively reliable given the large size of the dataset used.  
