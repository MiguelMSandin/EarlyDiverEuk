#----
#---- loading packages ----

library(ape)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(treeio)
library(tidytree)
library(EnvStats)

#----
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")
rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()
#---- 
#---- Set names of trees, files and output ---------------------------------------------------------

(files = list(
	trees=grep("clade_((Alveolata)|(Amoebozoa)|(Archaeplastida)|(Cryptista)|(Discoba)|(Haptista)|(Holozoa)|(Metamonada)|(Nucletmycea)|(Rhizaria)|(Stramenopila))*\\.tre$", dir(recursive=TRUE), value=TRUE), # The extracted clades
	abundances="../../resources/tipNames_OTUreads_norm.tsv", # The abundances table
	codes="../../resources/tipNamesIDs.tsv", # The coded tip names
	outDir="fractions",  # The output directory where all diversity fractions and plots are going to be exported
	out="plots/fractions/fractions", # A prefix for the main summary output
	minTips=300)) # A minimum number of tips to consider or not a supergroup

(files$root = unique(sub("\\/.*", "", files$trees)))
(files$replicates = unique(files$trees %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)))
(files$clades = unique(sub(".*\\/", "", files$trees)))

for(i in files$root){cat("  ", i, "\t", sum(grepl(i, files$trees)), "\n")}
for(i in files$replicates){cat("  ", i, "\t", sum(grepl(i, files$trees)), "\n")}
for(i in files$clades){cat("  ", i, "\t", sum(grepl(i, files$trees)), "\n")}
rm(i)

if(!dir.exists(dirname(files$out))){dir.create(dirname(files$out))}

#---- 
#---- Estimate species richness based on reads abundances ------------------------------------------

rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()

abun = fread(files$abundances)
tmp = fread(files$codes, header=FALSE)
all(tmp$V1 %in% abun$V1 & abun$V1 %in% tmp$V1)
abun = merge(abun, tmp, by="V1"); rm(tmp)
colnames(abun) = c("id", "abun", "code")

fileout = file(paste0(files$out, "_GoF.tsv"), "a")
write(paste("file", "clade", "tips",
            "truncated", "truncatedTips", 
            "truncatedGoFstatistic", "truncatedGoFpvalue", "truncatedQQslope", "truncatedQQadjR2",
            "ml2log2", "ml2log2Tips", 
            "ml2log2GoFstatistic", "ml2log2GoFpvalue", "ml2log2QQslope", "ml2log2QQadjR2",
            sep="\t"), 
      file=fileout)
dataOut = data.frame()

for(root in files$root){
  roots = grep(root, files$trees, value=TRUE)
  for(replicate in files$replicates){
    clades = grep(replicate, roots, value=TRUE)
    
    # Create directories
    diro = paste0(root, "/", replicate, "/", files$outDir)
    if(!dir.exists(diro)){dir.create(diro)}
    pdf(paste0(diro, "/", files$outDir, "_preston_fit.pdf"), width=11.69, height=8.27, paper='special')
    
    for(treen in clades){
      cat("\r  Working on tree ", grep(treen, files$trees), "/", length(files$trees), ": ", treen, "\t\t          ", sep="")
      clade = treen %>% sub(".*_", "", .) %>% sub("\\.tre", "", .)
      
      # read main tree
      tree = read.tree(treen)
      
      # fit Preston's models
      if(!all(tree$tip.label %in% abun$code)){cat("\n    Warning, not all tips are found in the codes provided...\n")}
      
      # Get abundances
      treeAbun = subset(abun, code %in% tree$tip.label)$abun
      
      # Fit Preston model
      preston = list(model_oc = prestonfit(treeAbun), 
                     model_ll = prestondistr(treeAbun))
      # Estimate number of expected species
      preston$model_oc_ext = veiledspec(preston$model_oc)
      preston$model_ll_ext = veiledspec(preston$model_ll)
      
      # Check goodness of fit
      gof = list()
      gof$oc_gof = gofTest(preston$model_oc$freq, preston$model_oc$fitted)
      gof$ll_gof = gofTest(preston$model_ll$freq, preston$model_ll$fitted)
      # gof$oc_r2r = lm(preston$model_oc$fitted ~ preston$model_oc$freq)
      # gof$ll_r2r = lm(preston$model_ll$fitted ~ preston$model_ll$freq)
      gof$oc_qq = lm(sort(preston$model_oc$fitted) ~ sort(preston$model_oc$freq))
      gof$ll_qq = lm(sort(preston$model_ll$fitted) ~ sort(preston$model_ll$freq))
      
      # Create a data.frame for plotting
      data = data.frame(x=1:length(preston$model_oc$freq),
                        octaves=2^(1:length(preston$model_oc$freq)),
                        freq=preston$model_oc$freq, 
                        truncated=preston$model_oc$fitted, 
                        mlLog2=preston$model_ll$fitted)
      
      preston$truncated = Ntip(tree)/preston$model_oc_ext[1]*100
      preston$mlLog = Ntip(tree)/preston$model_ll_ext[1]*100
      
      prestonPlot = ggplot(data)+
        geom_bar(aes(x=x, y=freq), stat="identity", alpha=0.2)+
        geom_line(aes(x=x, y=truncated), colour="orangered3")+
        geom_line(aes(x=x, y=mlLog2), colour="steelblue3")+
        theme_minimal()+
        labs(title=paste0(clade, ": ", Ntip(tree), " tips"),
             subtitle=paste0("This tree represents an estimated ", 
                             round(preston$truncated, 2), "% (Truncated fit; orange) or ", 
                             round(preston$mlLog, 2), "% (MaxLikelihood to Log2; blue) of the global diversity.\n",
                             "  Goodness-of-fit Truncated: KS distance: ", round(gof$oc_gof$statistic, 4), 
                             ";\tQ-Q slope: ", round(gof$oc_qq$coeff[[2]], 4), " (adj R²: ", round(summary(gof$oc_qq)$adj.r.squared, 4), ")",
                             "\n  Goodness-of-fit MLlog2:     KS distance: ", round(gof$ll_gof$statistic, 4),
                             ";\tQ-Q slope: ", round(gof$ll_qq$coeff[[2]], 4), " (adj R²: ", round(summary(gof$ll_qq)$adj.r.squared, 4), ")"),
             y="Species", x="Preston lognormal model frequencies observed")+
        scale_x_continuous(breaks=data$x, labels=data$octaves)
      plot(prestonPlot)
      
      lineout = paste(treen, clade, Ntip(tree), 
                      preston$truncated[[1]], preston$model_oc_ext[[1]], 
                      gof$oc_gof$statistic, gof$oc_gof$p.value, gof$oc_qq$coeff[[2]], summary(gof$oc_qq)$adj.r.squared,
                      preston$mlLog[[1]], preston$model_ll_ext[[1]], 
                      gof$ll_gof$statistic, gof$ll_gof$p.value, gof$ll_qq$coeff[[2]], summary(gof$ll_qq)$adj.r.squared, 
                      sep="\t")
      write(lineout, file=fileout, append=TRUE)
      dataOut = rbind(dataOut, data.frame(tree=treen,
                                          clade=clade,
                                          tips=Ntip(tree),
                                          truncated=preston$truncated[[1]],
                                          truncatedTips=preston$model_oc_ext[[1]],
                                          truncatedGoFstatistic=gof$oc_gof$statistic,
                                          truncatedGoFpvalue=gof$oc_gof$p.value,
                                          truncatedQQslope=gof$oc_qq$coeff[[2]],
                                          truncatedQQadjR2=summary(gof$oc_qq)$adj.r.squared,
                                          ml2log2=preston$mlLog[[1]],
                                          ml2log2Tips=preston$model_ll_ext[[1]],
                                          ml2log2GoFstatistic=gof$ll_gof$statistic,
                                          ml2log2GoFpvalue=gof$ll_gof$p.value,
                                          ml2log2QQslope=summary(gof$ll_qq)$adj.r.squared,
                                          ml2log2QQadjR2=gof$ll_qq$coeff[[2]]))
    }
    dev.off()
  }
}; rm(root, roots, replicate, clades, diro, treen, tree, treeAbun, preston, data, prestonPlot, lineout)

# write.table(dataOut, paste0(files$out, "_preston_GoF.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
#---- Extract single fractions file ----------------------------------------------------------------

data = fread(paste0(files$out, "_preston.tsv"))

out = data.frame()
for(group in unique(data$clade)){
	tmp = subset(data, clade==group)
	out = rbind(out, data.frame(clade=group, div=(max(tmp$frac)/100)))
	out = rbind(out, data.frame(clade=group, div=(min(tmp$frac)/100)))
}; rm(group, tmp)

# write.table(data, paste0(files$outDir, "/fractions.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
