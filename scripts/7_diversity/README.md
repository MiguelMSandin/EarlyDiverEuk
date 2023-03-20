# 7_diversity/

Here we will calculate the fraction of the diversity that our trees correspond to.  
The estimation of the total number of eukaryotic *species* is a debated topic of research. Such estimates goes from 1.8-2 million *species* ([Costello et al., 20212](https://academic.oup.com/sysbio/article/61/5/871/1734723)) to 5±3 million *species* ([Costello et al., 2013](https://www.science.org/doi/10.1126/science.1230318)) or with narrower estimates pointing at 8.7 ± 1.3 million *species* ([Mora et al., 2011](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001127)).  
Here we will use different and simpler approaches to provide a maximum and minimum boundaries while acknowledging the difficulty in arriving to such estimates. To do so, we will define *species* as 99% OTUs of the 18S rDNA gene sequences and therefore avoiding any further ambiguity in the concept of species (see [Queiroz (2007)](https://academic.oup.com/sysbio/article/56/6/879/1653163) or [Zachos (2016)](https://link.springer.com/book/10.1007/978-3-319-44966-1) for further details on the so-called *species* problem).  
  
[0_normalizeReads.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/7_diversity/0_normalizeReads.R):  
In order to compare the different datasets, we have to normalize the reads, and we will do so at 1000 reads. This arbitrary number was chosen on the basis of the OTU representative of most sequences from the smallest dataset. In our case, PR2 has one OTU representative of 1418 other sequences and PacBio has an OTU with 420101 reads. Therefore, 1000 is the closest power of 10 to the lowest maximum.
  
[1_diversityEstimate.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/7_diversity/1_diversityEstimate.R):  
In this script we will:  
- transform the scale of the trees to millions of years (if not has already been done manually in e.g. FigTree while visualizing them).  
- Export annotated supergroups from both the phylogenetic tree and the dated trees.  
- Estimate the total phylogenetic diversity based on the rarefied phylogenetic trees (taking into account the phylogenetic distance).  
- Estimate the total species richness based on the Preston's Lognormal Model ([Preston, 1948](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/1930989)) and the truncated lognormal model, implemented in the package '[vegan](https://github.com/vegandevs/vegan)'.
  
[2_globalPatterns.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/7_diversity/2_globalPatterns.R):  
Lastly, we will compare all different estimates of all supergroups across all the phylogenetic trees.  
  
Summarising, we found that our 99% OTU of the 18S rDNA trees represent in average a maximum of 74.8 (±11.8) % and an average minimum of 40.7 (±11.4) % of the *total* diversity. Amoebozoa showed the highest maximum with an average of 88.4 (±1.48) % and Cryptista showed the lowest minimum fraction at 13.4 (±0.15) %:  
  
|clade|min|max|
|---|---|---|
|Discoba|51.2 (±0.13) %|83.6 (±0.09) %|
|Metamonada|50.3 (±0.27) %|81.1 (±0.36) %|
|Amoebozoa|53.1 (±0.85) %|88.4 (±1.48) %|
|Nucletmycea|43.4 (±0.13) %|84.1 (±0.05) %|
|Holozoa|47.0 (±0.04) %|82.6 (±0.05) %|
|Haptista|49.8 (±0.06) %|79.6 (±0.03) %|
|Cryptista|13.4 (±0.15) %|44.1 (±0.06) %|
|Archaeplastida|41.5 (±0.14) %|69.9 (±0.23) %|
|Alveolata|31.5 (±0.13) %|68.5 (±0.09) %|
|Rhizaria|32.6 (±0.38) %|71.6 (±0.25) %|
|Stramenopila|34.3 (±0.19) %|70.0 (±0.10) %|

