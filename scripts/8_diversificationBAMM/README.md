# 8_diversificationBAMM

In this step we will try to identify shifts in diversification rates from the time-calibrated trees obtained from step [5_phyloDating/](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/5_phyloDating/) and the diversity fractions estimated in step [7_diversity](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/7_diversity/) with the [BAMM](10.1371/journal.pone.0089543) and [BAMMtools](10.1111/2041-210X.12199) packages.

```
```

[1_estimateBestShiftConfiguration.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/1_estimateBestShiftConfiguration.R):  
Firstly, we will run preliminary diversification analyses over 4 randomly selected trees to **estimate the most plausible number of expected shifts per supergroup** by arbitrarily setting the number of shifts to 10 and 50. To do so, we will generate the control file with the function '*generateControlFile*' and the priors estimated by the functions '*setBAMMpriors*' implemented in the [BAMMtools](10.1111/2041-210X.12199) R package. We will run the MCMC chain for 10 000 000 generations and setting the diversity fraction to the minimum and maximum estimates for every supergroup.  
  
[2_BAMMestimateBestShifts.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/2_BAMMestimateBestShifts.sh):  
Here we will run BAMM on the previously generated control files at the same time with the launcher script [2.1_BAMMestimateBestShifts_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/2.1_BAMMestimateBestShifts_launcher.sh). You might want to modify the headings of these scripts if you are running it in a cluster (see comments on the scripts).  
  
[3_getBestShiftConfiguration.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/3_getBestShiftConfiguration.R):  
Once all preliminary BAMM analyses have finished, we proceed to estimate the most plausible number of expected shifts based on the posterior results, or in BAMM terminology *the best shift configuration*.  
  
When estimating the best shift configuration we can also explore the effective sample size (ESS). And we see that the groups containing the largest number of tips (Alveolata, Holozoa and Nucletmycea) the **ESS is much lower than 200** (considered to be a threshold of convergence). Increasing the number of generations at the next script could increase EES, however downstream analyses will be exponentially affected and require very large amounts of RAM memory (>300 GB). Therefore we will cope with this issue by **using the relatively large number of phylogenetic trees as *replicates*** and discard shifts in diversification rates happening in a low number of *replicates*.  
  
[4_BAMMgenerateControl.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/4_BAMMgenerateControl.R):  
  
  
[5_BAMM.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/5_BAMM.sh):  
[5.1_BAMM_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/5.1_BAMM_launcher.sh).  
  
[6.2_rttBAMMextract.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6.2_rttBAMMextract.R):  
[6.1_BAMMextractRTT_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6.1_BAMMextractRTT_launcher.sh).  
  
[6_BAMMextractRTT.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6_BAMMextractRTT.sh):  
  
[7_analyzeBAMM.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/7_analyzeBAMM.R):  
  
[8_analyzeBAMM_all.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/8_analyzeBAMM_all.R):  
  




