#---- 
#---- loading packages ----

library(ape)
library(BAMMtools)
library(data.table)
library(dplyr)

#----
# 'setwd' to the root directory of a dated tree, where you have the subclades and the control files at the very least
setwd("")
#---- 
#---- BAMM set final control file for all clades, diversity estimates and proper Nshifts -----------

rm(list=ls()[!ls() %in% c("files")])

files = list(dirTrees="clades/",
			 divFile=paste0("fractions/", grep("fractions\\.tsv$", dir("fractions/"), value=TRUE)),
			 shiftsFile="output/shifts/shifts.tsv",
			 generations="10000000",
			 dirCtlOut="bamm/controls",
			 dirDivOut="bamm/diver")
files$trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE)
files

if(!dir.exists(files$dirCtlOut)){dir.create(files$dirCtlOut)}
if(!dir.exists(files$dirDivOut)){dir.create(files$dirDivOut)}

div = fread(files$divFile)
shi = fread(files$shiftsFile)

i= 0
for(tree in files$trees){
	clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
	cat("  Working on ", clade, " (", (i = i+1), "/", length(files$trees), ")\n", sep="")
	
	if(!dir.exists(paste0(files$dirCtlOut, "/", clade))){dir.create(paste0(files$dirCtlOut, "/", clade))}
	if(!dir.exists(paste0(files$dirDivOut, "/", clade))){
		dir.create(paste0(files$dirDivOut, "/", clade))
	}else{
		dir.create(paste0(files$dirDivOut, "/", clade, "/others"))
		system(paste0("mv ", files$dirDivOut, "/", clade, "/diversification_* ", files$dirDivOut, "/", clade, "/others/."))
		}
	
	phylo = read.tree(paste0(files$dirTrees, "/", tree))
	estimates = c(min(div$frac[grep(clade, div$clade)]), max(div$frac[grep(clade, div$clade)]))
	totalTips = c(min(div$total[grep(clade, div$clade)]), max(div$total[grep(clade, div$clade)]))
	shift = round(shi$mean[grep(as.character(clade), shi$clade)])
	
	for(j in 1:2){
		e = estimates[j]
		t = totalTips[j]
		ctl = paste0(files$dirCtlOut, "/", clade, "/control_diversification_", clade, "_div", round(e), "_shifts", shift,".ctl")
		cat("   Generating", ctl, "\n")
		priors <- setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
		generateControlFile(file=ctl, type = 'diversification',
							params = list(
								treefile = paste0("../../../", files$dirTrees, "/", tree),
								outName = paste0("diversification_", clade, "_div", round(e), "_shifts", shift),
								globalSamplingFraction = e/100,
								numberOfGenerations = files$generations,
								mcmcWriteFreq = '100',
								eventDataWriteFreq = '100',
								printFreq = '100',
								overwrite = '1',
								seed = '-1',
								deltaT = '0.01',
								lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
								lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
								muInitPrior = as.numeric(priors['muInitPrior']),
								# minCladeSizeForShift=Ntips(phylo)*0.1,
								expectedNumberOfShifts = shift))
	}
}; rm(i, tree, clade, phylo, estimates, totalTips, shift, j, e, t, ctl, priors)

#----  
