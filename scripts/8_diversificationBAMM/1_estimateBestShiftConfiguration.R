#---- 
#---- loading packages ----

library(ape)
library(BAMMtools)
library(data.table)
library(dplyr)

#----
setwd("")
# .rs.restartR()
#----
#---- Set file and directory names and select 4 random trees ---------------------------------------

# This script assumes that you have stored the trees in hierarchical directories
# For example: 'rootD/forward/RAcat1-1/step3f_RAcat1_rep1_rootD_renamed_treePLdated_median.tre' or 'rootD/forward/RAgamma1-2/step3f_RAgamma1_rep2_rootD_renamed_treePLdated_median.tr'
# And that in each directory you have the dated tree saved in newick format, a folder for the supergroup clades and another for the fraction files (as an output from the script '7_diversity/1_diversityEstimate.R')

files = list(trees = grep("rootD.*_median\\.newick", dir(recursive=TRUE), value=TRUE),
			 outDir="bamm",
			 outControl="bamm/control_test",
			 prefixControl="controlTest_diversification_",
			 prefixControlOut="diversification_",
			 dirTrees="clades",
			 pathTreesRel="../../../",
			 dirFractions="fractions",
			 shifts=c(10, 50),
			 generations="10000000")

trees = sample(files$trees, size=4)

#----
#---- BAMM set initial control file for all subclades and given diversity estimate -----------------

for(tree in trees){
	cat("  Working on tree ", tree, "\n")
	dir = sub("\\/[^\\/]+$", "", tree)
	if(!dir.exists(paste0(dir, "/", files$outDir))){dir.create(paste0(dir, "/", files$outDir))}
	if(!dir.exists(paste0(dir, "/", files$outControl))){dir.create(paste0(dir, "/", files$outControl))}
	
	div = grep("fractions\\.tsv$", dir(paste0(dir, "/", files$dirFractions)), value=TRUE)
	div = fread(paste0(dir, "/", files$dirFractions, "/", div))
	
	clades = grep("^clade.*.tre$", dir(paste0(dir, "/", files$dirTrees)), value=TRUE)
	
	i = 0
	for(clade in clades){
		group = clade %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
		cat("    Working on ", clade, " (", (i = i+1), "/", length(clades), ")\n", sep="")
		
		if(!dir.exists(paste0(dir, "/", files$outControl, "/", group))){dir.create(paste0(dir, "/", files$outControl, "/", group))}
		
		phylo = read.tree(paste0(dir, "/", files$dirTrees, "/", clade))
		
		estimates = c(min(div$frac[grep(group, div$clade)]), max(div$frac[grep(group, div$clade)]))
		totalTips = c(min(div$total[grep(group, div$clade)]), max(div$total[grep(group, div$clade)]))
		
		for(j in 1:2){
			e = estimates[j]
			t = totalTips[j]
			for(s in files$shifts){
				ctl = paste0(dir, "/", files$outControl, "/", group, "/", files$prefixControl, group, "_div", round(e), "_shifts", s,".ctl")
				cat("     Generating", ctl, "\n")
				priors <- setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
				generateControlFile(file=ctl, type = 'diversification',
									params = list(
										treefile = paste0(files$pathTreesRel, files$dirTrees, "/", clade),
										outName = paste0(files$prefixControlOut, group, "_div", round(e), "_shifts", s),
										globalSamplingFraction = e/100,
										numberOfGenerations = files$generations,
										mcmcWriteFreq = '100',
										eventDataWriteFreq = '100',
										printFreq = '100',
										overwrite = '1',
										seed = '-1',
										# simulatePriorShifts = '1',
										deltaT = '0.01',
										lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
										lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
										muInitPrior = as.numeric(priors['muInitPrior']),
										# minCladeSizeForShift=Ntips(phylo)*0.1,
										expectedNumberOfShifts = s))
			}
		}
	}
}; rm(tree, dir, clade, clades, group, i, phylo, estimates, totalTips, j, e, t, s, ctl, priors)

#----
