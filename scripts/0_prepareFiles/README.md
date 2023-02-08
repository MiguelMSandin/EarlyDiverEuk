# 0_prepareFiles/

Here we simply want to download and assemble the datasets for downstream analyses.  
  
Create a directory and '''cd''' to it. For example:  
'''bash
mkdir eukEcoEvo
cd eukEcoEvo
'''
Now create another directory for the 'data' and again '''cd''' to it:  
'''bash
mkdir data
cd data
'''

[0_download.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/0_prepareFiles/0_download.sh):  
With this script we will download all in-house dependencies and the original datasets. We will also modify the sequence names of the fasta files to accommodate different software requirements.  

[1_clusterPR2.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/0_prepareFiles/1_clusterPR2.sh):  
Once downloaded and before merging the two datasets, we will cluster the PR2 dataset into 99% OTUs to avoid redundancy.  

[2_mergeDataSets.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/0_prepareFiles/2_mergeDataSets.sh):  
Then, we remove those sequences of the PR2 dataset that are 100% identical to the PacBio dataset and finally we merge the two datasets.  

[3_subsetFinal.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/0_prepareFiles/3_subsetFinal.sh):  
And now we subset the whole dataset based on biological relevance. So we have 3 files: The complete dataset, a file gathering OTUs representative of at least 10 other sequences and a file gathering OTUs representative of at least 2 other sequences.
Since there are some critical groups with few sequences in databases, we will manually complement the two latter files with such sequences (from the lists [otus_reads10_complemented.list](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/resources/otus_reads10_complemented.list) and [otus_reads2_complemented.list](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/resources/otus_reads2_complemented.list)). These groups are: Ancoracystida, Cephalochordata, Filasterea, Glaucocystophyceae, Gromia, Hemichordata, Katablepharidaceae, Malawimonadidae, Mantamonadida, Mesostigmatophyceae, Monothalamids, Noctilucophyceae, Palmophyllophyceae, Parabasalia, Placozoa, Pluriformea, Preaxostyla, Protalveolata, Rhodelphea, Synchromophyceae and Tubothalamea.

Before beginning the analyses, it is important to keep a clear and structured directory. I preferred to let this aspect to your personal choice. Although I opted for the following structure:  
'''bash
cd eukEcoEvo
tree -L 2 -d
├── data
├── output
├── phylo
│   ├── step1f
│   ├── step1r
│   ├── step2f
│   ├── step2r
│   ├── step3f
│   └── step3r
├── phylo_dated
│   ├── rootD
│   └── rootA
├── plots
├── resources
└── scripts
    ├── 0_prepareFiles
    ├── 1_constraintTree
    ├── 2_phyloStep1
    ├── 3_phyloStep2
    ├── 4_phyloStep3
    ├── 5_phyloDating
    ├── 6_processDatedTrees
    ├── 7_diversity
    ├── 8_diversificationBAMM
    └── 9_diversificationCLaDS
'''
