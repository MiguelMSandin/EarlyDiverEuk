# 8_diversificationBAMM

In this step we will try to identify shifts in diversification rates using [BAMM](http://bamm-project.org/) and [BAMMtools](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12199) packages over the time-calibrated trees obtained from step [5_phyloDating/](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/5_phyloDating/) and the diversity fractions estimated in step [7_diversity](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/7_diversity/).
  
  
[1_estimateBestShiftConfiguration.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/1_estimateBestShiftConfiguration.R):  
Firstly, we will run preliminary diversification analyses over 4 randomly selected trees to **estimate the most plausible number of expected shifts per supergroup** by arbitrarily setting the number of shifts to 10 and 50. To do so, we will generate the control file with the function '*generateControlFile*' and the priors estimated by the functions '*setBAMMpriors*' implemented in the BAMMtools R package. We will run the MCMC chain for 10 000 000 generations and setting the diversity fraction to the minimum and maximum estimates for every supergroup.  
  
[2_BAMMestimateBestShifts.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/2_BAMMestimateBestShifts.sh):  
Here we will run BAMM on the previously generated control files at the same time with the launcher script [2.1_BAMMestimateBestShifts_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/2.1_BAMMestimateBestShifts_launcher.sh). You might want to modify the headings of these scripts if you are running it in a cluster (see comments on the scripts).  
  
[3_getBestShiftConfiguration.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/3_getBestShiftConfiguration.R):  
Once all preliminary BAMM analyses have finished, we proceed to estimate the most plausible number of expected shifts based on the posterior results from the `*mcmc_out.txt.gz` files with the functions `computeBayesFactors` and `plotPrior` from the R package BAMMtools.  
  
When estimating the best shift configuration we can also explore the [Effective Sample Size](https://beast.community/ess_tutorial) (ESS) with the function `effectiveSize` from the R package [coda](https://cran.r-project.org/web/packages/coda/index.html). And we see that the groups containing the largest number of tips (Alveolata, Holozoa and Nucletmycea) the **ESS is much lower than 200** (considered to be a threshold for convergence). Increasing the number of generations at the next script could increase EES, however downstream analyses will be exponentially affected and require very large amounts of RAM memory (>300 GB). Therefore we will cope with this issue by **using the relatively large number of phylogenetic trees as *replicates*** and discard shifts in diversification rates happening in a low number of *replicates*.  
  
[4_BAMMgenerateControl.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/4_BAMMgenerateControl.R):  
With this script we will generate the **final control files** to be run in BAMM. As in the script [1_estimateBestShiftConfiguration.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/1_estimateBestShiftConfiguration.R), we will generate the control file with the R package BAMMtools using the functions `generateControlFile` and `setBAMMpriors`. We will run the MCMC chain for 10 000 000 generations and setting the diversity fraction to the minimum and maximum estimates for every supergroup.  
  
[5_BAMM.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/5_BAMM.sh):  
[5.1_BAMM_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/5.1_BAMM_launcher.sh).  
This script will allow us to launch all **final BAMM analyses**, as done with the script [2_BAMMestimateBestShifts.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/2_BAMMestimateBestShifts.sh).  
  
[6_BAMMextractRTT.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6_BAMMextractRTT.sh):  
Once all BAMM analyses are complete, we will extract from every `*event_data.txt.gz` files the diversification Rates Through Time (RTT; for visualization purposes only) tables and **the best shift configuration** data with the functions `getEventData` and `getBestShiftConfiguration` from the R packages BAMMtools. To do so, we will loop through every `*event_data.txt.gz` file, and launch the script [6.1_BAMMextractRTT_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6.1_BAMMextractRTT_launcher.sh) which calls the R script [6.2_rttBAMMextract.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6.2_rttBAMMextract.R).  
This step is purposely convoluted to allow modifying headings according to different cluster needs and to ease the relaunch of possible failures.  
  
[7_analyzeBAMM.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/7_analyzeBAMM.R):  
Once we have the tables from the previous script, we can proceed to plotting them. With this script we plot the best shift configuration (and RTT slopes) of all supergroups from one time-calibrated tree of eukaryotes.  
  
[8_analyzeBAMM_all.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/8_analyzeBAMM_all.R):  
And if you want to explore the best shift configuration (and RTT slopes) of all time-calibrated trees of eukaryotes at once, this script will do so.  
  
