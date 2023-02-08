# 1_constraintTree/

In this step we will build the initial constraint tree:  
First, write manually (or download) a newick file that will be used as a backbone to include the sequence identifiers from the fasta file as tree tips. This first tree should look like:  
```
(A,((B,C),(D,(E,F))));
```
Where each letter is any given group (e.g.; Radiolaria, Diatomea, Dinoflagellata, ...).  
Then, follow the comments of the following script:  
  
[1_buildConstraintTree.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/1_constraintTree/1_buildConstraintTree.sh):  
Here you will check clade by clade that you are actually build a comprehensive constraint tree.  
This step is **the most laborious**. If you are building the tree from scratch, you might have to come back to this step after the first phylogenetic analyses in order to remove constraint that produce near-0 branch-lengths or polytomies.
