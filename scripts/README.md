# Scripts

In this folder you will find all scripts needed to replicate the analyses in a methodological order. You can find a more detailed readme file in each subfolder. Briefly:

## 0_prepareFiles/
Download and prepare files for downstream analyses.

## 1_constraintTree/
Build the initial constraint tree.

## 2_phyloStep1/
Build the initial phylogenetic tree including only OTUs with 10 or more reads using Maximum Likelihood in RAxML (GTR+CAT and GTR+Gamma) and RAxML-ng (GTR+Gamma).

## 3_phyloStep2/
Continue building the phylognetic tree including OTUs with 2 or more reads using maximum likelihood in IQTree, and the previous trees as constraint.

## 4_phyloStep3/
Finish building the phylogenetic tree including all OTUs using maximum likelihood in IQTree, and the previous trees as constraint.

## 5_phyloDating/
Calibrate in time the phylogenetic trees using Maximum Penalized likelihood to date the phylogenetic distance in TreePL.

## 6_processDatedTrees/
Extract Lineages Through Time plots and slopes of all dated trees and extract sub-clades of most abundant supergroups (with more than 300 tips).

## 7_diversity/
Estimate the fraction of the total diversity that each tree represent by different approaches (Preston Log-Normal model and an in-house approach based on phylogenetic distance).

## 8_diversificationBAMM/
Estimate diversification analyses in BAMM to extract potential shifts in diversification rates.

## 9_diversificationCLaDS/
Estimate diversification analyses in CLaDS to extract diversification rates.
