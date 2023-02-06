# Resources/

In this directory you will find several resources to help you replicate the analyses.
Some of which are created automatically along the pipeline, others are input files and others are log files to help keeping track of the different steps and outputs.

**EukProt_v9.tre**:
The original constraint tree used in this study, taken form the [EukProt](https://github.com/beaplab/EukProt) (and here the citation: [Richter *et al*., 2022](https://peercommunityjournal.org/articles/10.24072/pcjournal.173/)) database and modified to accommodate the different taxonomic framework used in [PR2]()https://pr2-database.org/ and the [PacBio](https://figshare.com/articles/dataset/Global_patterns_and_rates_of_habotat_transitions_across_the_eukaryotic_tree_of_life/15164772) sequences.

**info_trees.tsv**:
A table gathering relevant log information on the progressive phylogenetic reconstruction approach, such as the software used at every analysis, the constraint, the model, RAM memory used, CPU time, the likelihood of the resulting tree, number of tips, number of long branches pruned, number of intruders removed and number of tips after cleaning the raw tree.

**node_dates.tsv**:
The original and raw table compiled to choose the calibration nodes. This table contains all initially considered calibration nodes, their references and whether it has been used or not in this study.

**otu_number_lineage.tsv**:
A table gathering the number of tips per phylogenetic analysis and main lineage (e.g.; Discoba, CRuMs, Hemimastigophora, Rhizaria, Nucletmycea).

**otu_reads.tsv**:
A two columns table with the sequence identifier and the number of reads.

**otus_reads10_complemented.list**:
A list (or one column table) with all sequence identifiers from OTUs that have 10 or more reads.

**otus_reads2_complemented.list**:
A list with all sequence identifiers from OTUs that have 2 or more reads.

**tipNamesIDs.tsv**:
A two columns table with the sequence identifier and a short identifier to reduce the size of the trees.

**tipNamesIDs_reverse.tsv**:
The reverse of the previous table (short identifier and the original sequence identifier) to add again the original names to the reduced trees.
