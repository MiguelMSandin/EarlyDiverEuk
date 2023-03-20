# Eco-evolution of eukaryotes
In this repository you will find all materials, resources and scripts needed to replicate the study:
***INSERT_LINK_HERE***

## Structure of the repository
In every directory of this repository, you will find a readme file explaining the contents of the given directory.
Briefly:
- [Resources/](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/resources): Contains both input files needed to carry out the analyses and log files from the original study.
- [Scripts/](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts): all scripts used in this study ordered methodologically

## Dependencies
- [BAMM](http://bamm-project.org/)
- [ClaDS](https://hmorlon.github.io/PANDA.jl/dev/)
- [Julia](https://julialang.org/):
  - Packages: ArgParse, JLD2, PANDA
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [mothur](https://mothur.org/)
- [Python](https://www.python.org/):
  - Modules: argparse, Bio (SeqIO, Phylo), ete3 (Tree), numpy, or, re, statistics, subprocess, sys
- [R](https://www.r-project.org/):
  - Packages: ape, BAMMtools, coda, data.table, dplyr, ggplot2, ggtree, HDInterval, optparse, phangorn, RPANDA, tidyr, treeio, vegan
- [RAxML](https://github.com/stamatak/standard-RAxML)
- [RAxML-ng](https://github.com/amkozlov/raxml-ng)
- [treeannotator](https://beast.community/treeannotator)
- [TreePL](https://github.com/blackrim/treePL)
- [trimAl](http://trimal.cgenomics.org/downloads)
- [IQ-TREE](http://www.iqtree.org/)
#### Optional
- [Rstudio](https://rstudio.com/products/rstudio/download/)
- [figTree](http://tree.bio.ed.ac.uk/software/figtree/)

### In-house dependencies (not included in this repository)
- [buildConstrainTree.sh](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/buildConstrainTree.sh)
- [checkConstrainTaxa.sh](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/checkConstrainTaxa.sh)
- [checkConstrainTree.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/checkConstrainTree.py)
- [fastaConcat.py](https://github.com/MiguelMSandin/random/blob/main/fasta/fastaConcat.py)
- [fastaRevCom.py](https://github.com/MiguelMSandin/random/blob/main/fasta/fastaRevCom.py)
- [fileNameReplace.py](https://github.com/MiguelMSandin/random/blob/main/others/fileNameReplace.py)
- [findSeqs.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/findSeqs.py)
- [newick2nexus.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/newick2nexus.py)
- [sequenceSelect.py](https://github.com/MiguelMSandin/random/blob/main/fasta/sequenceSelect.py)
- [treeCheckIntruders.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeCheckIntruders.py)
- [treeColourBranches.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeColourBranches.py)
- [treeCompareTips.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeCompareTips.py)
- [treeLTTsubTrees.R](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeLTTsubTrees.R)
- [treePruneList.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treePruneList.py)
- [treePruneOutliers.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treePruneOutliers.py)
- [treeRemoveBranchLengths.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeRemoveBranchLengths.py)
- [treeRootOutgroup.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeRootOutgroup.py)
- [treeTipRename.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeTipRename.py)
#### Optional
- [treeLCAcount.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeLCAcount.py)
- [treeStats.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/treeStats.py)
