# 6_processDatedTrees

In this step we will simply process the dated trees, and extract useful direct information from the node heights, such as the Lineages Through Time (LTT) slopes and the node ages for all major eukaryotic supergroups.

[1_processDatedTrees.sh](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/6_processDatedTrees/1_processDatedTrees.sh):  
With this script we will extract the LTT table and plots for each tree (with the ```Rscript treeLTTsubTrees.R```), rename the tree to the old and long human-readable tip names for manual checks and finally colour the branches for a quick visual inspection.  
  
[2_getLTTsummary.R](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/6_processDatedTrees/2_getLTTsummary.R):  
This script will combine and summarize all LTT tables. It is recommended to open it in Rstudio (or alike), so you can modify parameters and aesthetics on the go.  
  
[3_getAges.R](https://github.com/MiguelMSandin/EukEcoEvo/blob/main/scripts/6_processDatedTrees/3_getAges.R):  
Lastly, this script will loop through all trees and extract the appearance age of all major eukaryotic supergroups for a quick comparison among phylogenetic trees and the two different rooting scenarios (in Discoba and Amorphea).  
  
