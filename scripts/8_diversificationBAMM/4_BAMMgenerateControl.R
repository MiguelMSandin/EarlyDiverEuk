#---- 
#---- loading packages ----

library(ape)
library(BAMMtools)
library(data.table)
library(dplyr)

#----
# 'setwd' to the root directory of a dated tree, where you have the subclades and the control files at the very least
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC02-2/rootD/")
#---- 
#---- BAMM set final control file for all clades, diversity estimates and proper Nshifts -----------

folders = grep("RA", dir(), value=TRUE)
rm(list=ls()[!ls() %in% c("folders")])

for(folder in folders){
  cat("Setting working directory in ", folder, "\t", grep(folder, folders), "/", length(folders), "          \n", sep="")
  setwd(folder)
  
  files = list(dirTrees="clades/",
               divFile=grep("fractions\\.tsv$", dir(), value=TRUE),
               shiftsFile="../../../../output/shifts/shifts.tsv",
               generations="10000000",
               writeFreq="100",
               dirOut="bamm")
  files$trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE)
  files
  
  if(!dir.exists(files$dirOut)){dir.create(files$dirOut)}
  
  div = fread(files$divFile)
  shi = fread(files$shiftsFile)
  
  i= 0
  for(tree in files$trees){
    clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
    cat("  Working on ", clade, " (", (i = i+1), "/", length(files$trees), ")\n", sep="")
    
    if(!dir.exists(paste0(files$dirOut, "/", clade))){dir.create(paste0(files$dirOut, "/", clade))}
    
    phylo = read.tree(paste0(files$dirTrees, "/", tree))
    estimates = c(min(div$frac[grep(clade, div$clade)]), max(div$frac[grep(clade, div$clade)]))
    totalTips = c(min(div$total[grep(clade, div$clade)]), max(div$total[grep(clade, div$clade)]))
    shift = round(shi$mean[grep(as.character(clade), shi$clade)])
    factor = ifelse(clade =="Alveolata" | clade=="Holozoa" | clade=="Nucletmycea", 10, 1)
    
    for(j in 1:2){
      e = estimates[j]
      t = totalTips[j]
      ctl = paste0(files$dirOut, "/", clade, "/control_diversification_", clade, "_div", round(e), "_shifts", shift,".ctl")
      cat("   Generating", ctl, "\n")
      priors = setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
      generateControlFile(file=ctl, type = 'diversification',
                          params = list(
                            treefile = paste0("../../", files$dirTrees, "/", tree),
                            outName = paste0("diversification_", clade, "_div", round(e), "_shifts", shift),
                            globalSamplingFraction = e/100,
                            numberOfGenerations = format(as.numeric(files$generations) * factor, scientific=FALSE),
                            mcmcWriteFreq = format(as.numeric(files$writeFreq) * factor, scientific=FALSE),
                            eventDataWriteFreq = format(as.numeric(files$writeFreq) * factor, scientific=FALSE),
                            printFreq = format(as.numeric(files$writeFreq) * factor, scientific=FALSE),
                            overwrite = '1',
                            seed = '-1',
                            deltaT = '0.01',
                            lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
                            lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
                            muInitPrior = as.numeric(priors['muInitPrior']),
                            # minCladeSizeForShift=Ntips(phylo)*0.1,
                            expectedNumberOfShifts = shift))
    }
  }; rm(i, tree, clade, phylo, estimates, totalTips, shift, j, e, t, ctl, priors, factor)
  setwd("..")
  cat("\n")
}

#----  
