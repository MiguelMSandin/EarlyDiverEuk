# 9_diversificationCLaDS//
  
In this step we will apply the Data Augmentation method to estimate the branch specific speciation rates, the diversification Rates Through Time (RTT) and the augmented Diversity Through Time (DTT) slopes of time-calibrated phylogenetic trees implemented in the recently developed [ClaDS](https://hmorlon.github.io/PANDA.jl/dev/) package.  
  
[0_getFractions.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/0_getFractions.sh):  
Firstly we will extract the corresponding diversity fractions of every supergroup and create two tables for the minimum and maximum estimates.  
  
[1_CLaDS.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/1_CLaDS.sh):  
Here we will launch the ClaDS implementation over all supergroups from all eukaryotic time-calibrated trees for both the minimum and maximum diversity fractions. To do so, we will call the script [1.1_CLaDS_launcher.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/1.1_CLaDS_launcher.sh) which will call the actual ClaDS function with the script [1.2_clads.jl](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/1.2_clads.jl). This step will export an `Rdata` object that we will read and analyze in the next script.  
As for the script [6_BAMMextractRTT.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/8_diversificationBAMM/6_BAMMextractRTT.sh) from the BAMM analyses, this step is purposely convoluted to allow modifying headings according to different cluster needs and to ease the relaunch of possible failures.  
  
[2_RTT-CLaDSextract.sh](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/2_RTT-CLaDSextract.sh):  
Once ClaDS have been run in all trees of interest, we will extract the RTT slopes, the DTT slopes and the branch specific speciation rates. We will loop through `*Rdata` file and call the script [2.1_rttClaDSextract.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/2.1_rttClaDSextract.R).
  
[3_analyzeCLaDSoutput.R](https://github.com/MiguelMSandin/EukEcoEvo/tree/main/scripts/9_diversificationCLaDS/3_analyzeCLaDSoutput.R):  
Lastly, we will analyze all different RTT and DTT slopes in an attempt to find patterns in the macro-diversification of eukaryotes.  
  
