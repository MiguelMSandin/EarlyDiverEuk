#---- 
#---- loading packages ----

library(ape)
library(treeio)

library(BAMMtools)
library(HDInterval)

library(data.table)
library(tidyr)
# library(dplyr)

library(ggplot2)
library(ggtree)

# library(devtools)
# install_github("hmorlon/PANDA", dependencies = TRUE)
library(RPANDA)

library(vegan)
library(drc)
library(nlme)
library(aomisc)

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, 
                          -1600 , -1400, -1200, 
                          -1000, -720, -635, 
                          -541, -485.4, -443.8, -419, -358.9, -298.9, 
                          -251.9, -201.4, -145, 
                          -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", 
                           "Calymnian","Ectasian","Stenian",
                           "Tonian","Cryogenian","Ediacran",
                           "Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian",
                           "Triassic","Jurassic","Cretaceous",
                           "Paleogene","Neogene","Quaternary"),
                  era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic",
                        "Mesoproterozoic","Mesoproterozoic","Mesoproterozoic",
                        "Neoproterozoic","Neoproterozoic","Neoproterozoic",
                        "Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
                        "Mesozoic","Mesozoic","Mesozoic",
                        "Cenozoic","Cenozoic","Cenozoic"))

#----  
# setwd("/home/miguel/Desktop/PhD/0_Thesis/2_chapter/data/diver/MC7_mcmctree/")
setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/RAcat1-1/")
setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/RAgamma1-1/")
rm(list=ls()[!ls() %in% c("geo")])
# .rs.restartR()
#---- 
#----  Transform scale of tree if needed and export independent clades -----------------------------

files = list(tree="step3f_RAgamma1_rep1_rootD_renamed_treePLdated_median.tre",
             factor=100,
             dirClades="clades",
             minTips=300)

# Create directory if doesn't exist
if(!dir.exists(files$dirClades)){dir.create(files$dirClades)}else{cat("  directory '", files$dirClades, "' already exists\n", sep="")}

# Read annotated nexus tree file
treei <- read.beast(files$tree)

hist(log(treed$branch.length))
summary(treed$branch.length)

# Transform scale if needed
# tree = read.tree(file)
# tree$edge.length = tree$edge.length * factor
# write.tree(tree, sub("\\.[^\\.]+$", "_timeScaled.tre", file))

# Transform tree into a data.frame
treed <- as_tibble(treei)
names(treed) <- gsub("!", "", names(treed))

# Select all annotated nodes
(annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))

# loop through all annotated clades, and export those with at least 'files$minTips' tips
for(annot in annotations){
    cat("  Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")\n", sep="")
    node = as.numeric(subset(treed, name==annot)$node)
    tmp = tree_subset(treei, node=node, levels_back=0)
    tmp = as.phylo(tmp)
    tmp$root.edge = NULL
    if(Ntip(tmp) >= files$minTips){
        cat("    Exporting tree with", Ntip(tmp),"tips\n")
        write.tree(tmp, paste0(files$dirClades, "/clade_", annot, ".tre"))
    }else{cat("    ! Clade not exported (", Ntip(tmp), " tips)\n", sep="")}
};rm(annot, node, tmp)

#

#----  Export independent clades from phylogenetic tree --------------------------------------------

files = list(tree="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/step3f/raxml-gamma/rootD/step3f_RAgamma1_iqtreef_GTRg_rep1_cleaned_rootD_coloured.tre",
             dirClades="clades/phylo",
             minTips=300)

# Create directory if doesn't exist
if(!dir.exists(files$dirClades)){dir.create(files$dirClades)}else{cat("  directory '", files$dirClades, "' already exists\n", sep="")}

# Read annotated nexus tree file
treei <- read.beast(files$tree)

# Transform tree into a data.frame
treed <- as_tibble(treei)
names(treed) <- gsub("!", "", names(treed))

# Select all annotated nodes
(annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))

# loop through all annotated clades, and export those with at least 'files$minTips' tips
for(annot in annotations){
    cat("  Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")\n", sep="")
    node = as.numeric(subset(treed, name==annot)$node)
    tmp = tree_subset(treei, node=node, levels_back=0)
    tmp = as.phylo(tmp)
    if(Ntip(tmp) >= files$minTips){
        cat("    Exporting tree with", Ntip(tmp),"tips\n")
        write.tree(tmp, paste0(files$dirClades, "/clade_", annot, ".tre"))
    }else{cat("    ! Clade not exported (", Ntip(tmp), " tips)\n", sep="")}
};rm(annot, node, tmp)

#

#----
#----
#---- Estimate phylogenetic species richness based on the rarefied trees raw -----------------------

rm(list=ls()[!ls() %in% c("geo")]); gc()
# .rs.restartR()

files = list(abundances="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNames_OTUreads_norm.tsv",
             useAbundance=TRUE,
             dirTrees="clades/phylo",
             outSlopes="RAgamma1-1_supergroups_rarefiedPhylo_fraction.pdf",
             outData=  "RAgamma1-1_supergroups_rarefiedPhylo_fraction.tsv",
             outDir="fractions",
             maxSamplingFactor=10) # X times the number of tips

# Create directory if doesn't exist
if(!dir.exists(files$outDir)){dir.create(files$outDir)}else{cat("  directory '", files$outDir, "' already exists\n", sep="")}

# Read abundances file
if(files$useAbundance){abun = fread(files$abundances)}

# Rarefy the trees _________________________________________________________________________________
pdf(paste0(files$outDir, "/", files$outSlopes), width=11.69, height=8.27, paper='special')
outData = data.frame(); outSlopes = data.frame()
for(tree in grep("^clade_.*\\.tre$", dir(files$dirTrees), value=TRUE)){
    # tree="clade_Cryptista.tre"
    # tree="clade_Metamonada.tre"
    cat(" ", tree)
    clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
    phylo = read.tree(paste0(files$dirTrees, "/", tree))
    # plot(phylo, show.tip.label = FALSE)
    
    # Get distances
    dists = data.frame(tips=phylo$tip.label,
                       dists=node.depth.edgelength(phylo)[0:Ntip(phylo)])
    total = sum(dists$dists)
    
    # Get abundances and merge with distances
    abuns = abun[abun$V1 %in% phylo$tip.label]
    dists = merge(dists, abuns, by.x="tips", by.y="V1")
    dists$abunRel = dists$V2 / max(dists$V2)
    
    # Rarefy
    rare = data.frame()
    k = 0
    samples = round(seq(1, Ntip(phylo)*files$maxSamplingFactor, length.out = 100))
    for(i in samples){
        cat("\r  ", tree, "\t", (k = k + 1), "%", sep="")
        tmp = c()
        for(j in 1:100){
            if(files$useAbundance){
                tips = sample(dists$tips, i, replace=TRUE, prob=dists$abunRel)
            }else{
                tips = sample(dists$tips, i, replace=TRUE)
            }
            # maybe try by extracting tips with a probability vector of their distance to the root?
            # Or include abundance instead, as follows:
            tips = tips[!duplicated(tips)]
            tmp = c(tmp, sum(dists[dists$tips %in% tips,]$dists))
        }
        rare = rbind(rare, tmp)
    }; cat("\n")
    
    # Extract stats
    rares = data.frame(clade=clade,
                       sample=samples,
                       mean=apply(rare, 1, mean),
                       min=apply(rare, 1, min),
                       max=apply(rare, 1, max))
    outSlopes <- rbind(outSlopes, rares)
    
    # Estimate diversity
    div1 = subset(rares, sample >= Ntip(phylo))[1,]$mean
    frac = div1/total*100
    divt = (Ntip(phylo)*100)/frac
    
    # Apply an Asymptotic regression model
    # Y=a−(a−b)exp(−cX)
    # a : maximum attainable Y
    # b : Y at x=0
    # c : proportional to the relative rate of Y increase while X increases
    # model <- drm(rares$mean ~ rares$sample, fct = DRC.asymReg())
    # plateau = summary(model)$coefficients[3]
    # rares$fitted = fitted(model)
    # div = plateau - (plateau - summary(model)$coefficients[1]) * exp(-summary(model)$coefficients[2] * Ntip(phylo))
    # divtm = (Ntip(phylo)*100)/div
    
    outData = rbind(outData, data.frame(clade= clade, 
                                        tips=Ntip(phylo), 
                                        fraction=frac,
                                        total=divt))
    
    # plot
    slope = ggplot(rares, aes(x=sample, y=mean))+
        geom_line(aes(y=min), colour="lightgrey")+
        geom_line(aes(y=max), colour="lightgrey")+
        geom_hline(yintercept=total, color="springgreen3")+
        # geom_hline(yintercept=plateau, color="orangered3")+
        # geom_segment(aes(x=Ntip(phylo), y=0, 
        #                  xend=Ntip(phylo), yend=div), colour="orangered3")+
        # geom_segment(aes(x=0, y=div, 
        #                  xend=Ntip(phylo), yend=div), colour="orangered3")+
        geom_segment(aes(x=Ntip(phylo), y=0, 
                         xend=Ntip(phylo), yend=div1), colour="springgreen3")+
        geom_segment(aes(x=0, y=div1, 
                         xend=Ntip(phylo), yend=div1), colour="springgreen3")+
        # geom_line(aes(y=fitted), color="orangered3")+
        geom_line()+
        theme_minimal()+
        labs(title=paste0(clade, ": ", Ntip(phylo), " tips"),
             subtitle=paste0("This tree represents an estimated ~", round(frac), "% of the global diversity (",
                             "estimated to be ~", round(divt), " total tips) based on the total phylogenetic diversity.\n"),
                             # "Or ", round(div, 2), "% of the global diversity (",
                             # "estimated to be ", round(divtm), " total tips) based on the fitted asymptotic regression model.\n"),
             y="Mean phylogenetic distance to Root", x="Sample size")
    plot(slope)
}; rm(tree, phylo, rare, rares, dists, total, abuns, div1, frac, divt, clade, k, samples, i, j, tmp, slope); dev.off()
# rm(model, plateau, div, divtm);
write.table(outData, paste0(files$outDir, "/", files$outData), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(outSlopes, paste0(files$outDir, "/", gsub("\\.tsv", "_raw.tsv", files$outData)), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#

#---- Estimate species richness based on reads abundances ------------------------------------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")]); gc()

files = list(abundances="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNames_OTUreads_norm_IDs.tsv",
             dirTrees="clades",
             dirRare="preston",
             outSlopes="RAgamma1-1_supergroups_preston_fraction.pdf",
             outData=  "RAgamma1-1_supergroups_preston_fraction.tsv",
             outDir="fractions")

# Create directory if doesn't exist
if(!dir.exists(files$outDir)){dir.create(files$outDir)}else{cat("  directory '", files$outDir, "' already exists\n", sep="")}

abundances = fread(files$abundances)

trees = grep("^clade_.*\\.tre$", dir(paste0(files$dirTrees)), value=TRUE); i = 0
dataOut = data.frame()
pdf(paste0(files$outDir, "/", files$outSlopes), width=11.69, height=8.27, paper='special')
for(tree in trees){
    cat("  Working on ", tree, " (", (i = i + 1), "/", length(trees), ")\n", sep="")
    
    clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
    
    phylo = read.tree(paste0(files$dirTrees, "/", tree))
    abun = abundances[abundances$V1 %in% phylo$tip.label,]$V2
    
    model_oc = prestonfit(abun)
    model_ll = prestondistr(abun)
    
    data = data.frame(x=1:length(model_oc$freq),
                      octaves=2^(1:length(model_oc$freq)),
                      freq=model_oc$freq, 
                      qPois=model_oc$fitted, 
                      mlLog2=model_ll$fitted)
    
    model_oc_ext = veiledspec(model_oc)
    model_ll_ext = veiledspec(model_ll)
    
    qPois = Ntip(phylo)/model_oc_ext[1]*100
    mlLog = Ntip(phylo)/model_ll_ext[1]*100
    
    prestonPlot = ggplot(data)+
        geom_bar(aes(x=x, y=freq), stat="identity", alpha=0.2)+
        geom_line(aes(x=x, y=qPois), colour="orangered3")+
        geom_line(aes(x=x, y=mlLog2), colour="steelblue3")+
        theme_minimal()+
        labs(title=paste0(clade, ": ", Ntip(phylo), " tips"),
             subtitle=paste0("This tree represents an estimated ", 
                             round(qPois, 2), "% (Quasi-Poisson fit; orange) or ", 
                             round(mlLog, 2), "% (MaxLikelihood to Log2; blue) of the global diversity."),
             y="Species", x="Preston lognormal model frequencies observed")+
        scale_x_continuous(breaks=data$x, labels=data$octaves)
    
    plot(prestonPlot)
    
    dataOut = rbind(dataOut, data.frame(clade=clade,
                                        tips=Ntip(phylo),
                                        quasiPoisson=qPois,
                                        quasiPoissonTips=model_oc_ext[1],
                                        maxLikelihood2log2=mlLog,
                                        maxLikelihood2log2Tips=model_ll_ext[1]))
    
}; rm(trees, i, tree, clade, phylo, abun, model_oc, model_ll, model_oc_ext, model_ll_ext, qPois, mlLog, data, prestonPlot);dev.off()
write.table(dataOut, paste0(files$outDir, "/", files$outData), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




#---- Check and compare the different estimates ----------------------------------------------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

files = list(dir="fractions/",
             preston="RAgamma1-1_supergroups_preston_fraction.tsv",
             phylo="RAgamma1-1_supergroups_rarefiedPhylo_fraction.tsv",
             outData="RAgamma1-1_supergroups_fractions.tsv",
             outPlot="RAgamma1-1_supergroups_fractions.pdf")

tmp1 = fread(paste0(files$dir, "/", files$preston))
tmp2 = fread(paste0(files$dir, "/", files$phylo))

if(all(tmp1$clade == tmp2$clade) & all(tmp1$tips == tmp2$tips)){
    data = data.frame(clade=rep(tmp1$clade, 3),
                      tips=rep(tmp1$tips, 3),
                      method=c(rep("quasiPoisson", nrow(tmp1)), rep("maxLikelihood2log2", nrow(tmp1)), rep("rarefied", nrow(tmp1))),
                      frac=c(tmp1$quasiPoisson, tmp1$maxLikelihood2log2, tmp2$fraction),
                      total=c(tmp1$quasiPoissonTips, tmp1$maxLikelihood2log2Tips, tmp2$total))
    rm(tmp1, tmp2)
}else{cat("\nPlease check that the clades and tips in both files are in matching order\n\n")}

write.table(data, paste0(files$dir, "/", files$outData), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

unique(data$clade)
data$clade = factor(data$clade, levels=c("Discoba", "Metamonada",
                                         "Amoebozoa", "Nucletmycea", "Holozoa",
                                         "Haptista",  "Cryptista", "Archaeplastida",
                                         "Rhizaria", "Stramenopila","Alveolata"))

(summ = ggplot(data, aes(x=clade, y=frac))+
    geom_point(aes(colour=method, size=tips))+
    theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Summary of diversity estimates",
         y="Diversity fraction estimate (percentage of total)", x="Clades"))

pdf(paste0(files$dir, "/", files$outPlot), width=11.69, height=8.27, paper='special'); plot(summ); dev.off()



#----
##### Analyse BAMM output ##########################################################################
#----
#---- BAMM set control file for the main tree accounting for the given diversity estimate ----------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(tree="step3f_RAcat1_rep1_rootD_renamed_treePLdated_median_newick.tre",
             divFile="RAcat1-1_supergroups_fractions.tsv",
             shifts=c(1, 5, 10, 50),
             dirCtlOut="controls")

if(!dir.exists(files$dirCtlOut)){dir.create(files$dirCtlOut)}

div = fread(files$divFile)
divs = div %>% dplyr::group_by(method) %>% dplyr::summarise(mean=mean(frac))


phylo = read.tree(files$tree)

estimates = c(min(divs$mean), max(divs$mean))
totalTips = c(round(Ntip(phylo)*100/min(divs$mean)), round(Ntip(phylo)*100/max(divs$mean)))

for(j in 1:2){
    e = estimates[j]
    t = totalTips[j]
    for(s in files$shifts){
        ctl = paste0(files$dirCtlOut, "/control_diversification_all_div", round(e), "_shifts", s,".ctl")
        cat("   Generating", ctl, "\n")
        priors <- setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
        generateControlFile(file=ctl, type = 'diversification',
                            params = list(
                                treefile = paste0("../", files$tree),
                                outName = paste0("diversification_all_div", round(e), "_shifts", s),
                                globalSamplingFraction = e/100,
                                numberOfGenerations = '1000000',
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
        # system(paste("echo '\nsimulatePriorShifts = 1' >>", ctl))
    }
    
}
; rm(trees, i, tree, clade, phylo, estimates, totalTips, j, e, t, s, ctl, priors)


#---- BAMM set initial control file for all subclades and given diversity estimate -----------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirTrees="clades/",
             divFile="fractions/RAgamma1-1_supergroups_fractions.tsv",
             shifts=c(10, 50),
             generations="10000000",
             dirCtlOut="bamm/controls")

if(!dir.exists(sub("\\/.*", "", files$dirCtlOut))){dir.create(sub("\\/.*", "", files$dirCtlOut))}
if(!dir.exists(files$dirCtlOut)){dir.create(files$dirCtlOut)}

div = fread(files$divFile)

trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE); i= 0
for(tree in trees){
    clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
    cat("  Working on ", clade, " (", (i = i+1), "/", length(trees), ")\n", sep="")
    
    if(!dir.exists(paste0(files$dirCtlOut, "/", clade))){dir.create(paste0(files$dirCtlOut, "/", clade))}
    
    phylo = read.tree(paste0(files$dirTrees, "/", tree))
    
    estimates = c(min(div$frac[grep(clade, div$clade)]), max(div$frac[grep(clade, div$clade)]))
    totalTips = c(min(div$total[grep(clade, div$clade)]), max(div$total[grep(clade, div$clade)]))
    
    for(j in 1:2){
        e = estimates[j]
        t = totalTips[j]
        for(s in files$shifts){
            ctl = paste0(files$dirCtlOut, "/", clade, "/control_diversification_", clade, "_div", round(e), "_shifts", s,".ctl")
            cat("   Generating", ctl, "\n")
            priors <- setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
            generateControlFile(file=ctl, type = 'diversification',
                                params = list(
                                    treefile = paste0("../../../", files$dirTrees, "/", tree),
                                    outName = paste0("diversification_", clade, "_div", round(e), "_shifts", s),
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
            # system(paste("echo '\nsimulatePriorShifts = 1' >>", ctl))
        }
        
    }
}; rm(trees, i, tree, clade, phylo, estimates, totalTips, j, e, t, s, ctl, priors)


#----  
##### Analyze diversification shifts from BAMM MCMC output chain -----------------------------------

# Read file
library(coda)
files = list(mcmcDir="bamm/Rhizaria/",
             burn=25)

files$clade = files$mcmcDir %>% sub("\\/$", "", .) %>% sub(".*\\/", "", .)
files$outPlot = paste0("bamm/shifts_", files$clade, ".pdf")
files$outTable = paste0("bamm/shifts_", files$clade, ".tsv")
files$files = grep("mcmc", dir(files$mcmcDir, recursive=TRUE), value=TRUE)

# Start the loop
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
outTable = data.frame()
for(file in files$files){
    cat("  Plotting ", file, "\t(", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    shifts = as.numeric(file %>% sub(".*shifts", "", .) %>% sub("_.*", "", .))
    div = file %>% sub(".*div", "", .) %>% sub("_.*", "", .)
    
    # Read the MCMC file
    mcmcout = read.csv(paste0(files$mcmcDir, "/", file), header=T)
    # Remove burnin
    postburn = mcmcout[((files$burn/100)*nrow(mcmcout)):nrow(mcmcout), ]
    
    # Start a plot of 4 grids
    # dev.off()
    par(mar=c(2.5,2.5,1,1))
    layout(matrix(c(1,6,2,4, 1,6,3,5), ncol=2),
           heights=c(1, 1, 3, 3))
    plot.new()
    text(0.5,0.5, file, cex=2, font=2)
    
    # Check Effective Sampling Site
    essS = effectiveSize(postburn$N_shifts)
    essL = effectiveSize(postburn$logLik)
    
    # Check for convergence
    plot(mcmcout$logLik ~ mcmcout$generation, 
         main="Convergence", xlab="Generation", ylab="Log-Likelihood")
    
    # Plot the number of shifts 
    table(postburn$N_shifts) / nrow(postburn)
    hist(postburn$N_shifts,
         xlim=c(0, max(postburn$N_shifts)),
         main="Histogram of number of shifts occurrence", xlab="Number of shifts", ylab="Frequency")
    shiftsM = names(which.max(table(postburn$N_shifts)))
    
    # Compute Bayes Factor
    bayes = computeBayesFactors(mcmcout, expectedNumberOfShifts=shifts, burnin=files$burn/100)[,1]
    barplot(bayes, 
            main="Bayes Factor", xlab="Number of shifts", ylab="Bayes Factor")
    bayesM = names(which.max(bayes))
    
    # Plot Prior Vs Posterior
    plotPrior(postburn, expectedNumberOfShifts=shifts)
    
    # Add subtitle with the summary of the statistics
    plot.new()
    text(0.5, 0.5, paste0("ESS shifts: ", round(essS, 1) , "; ESS LogLik:", round(essL, 1), 
                          " - Most freq shift: ", shiftsM,
                          "; Shift with highest Bayes factor: ", bayesM), 
         cex=1.5, font=1)
    # Add row to table to export
    outTable = rbind(outTable, data.frame(file=file,
                                          diversity=div,
                                          shiftPrior=shifts,
                                          mostFreqShift=shiftsM,
                                          bayesShift=bayesM,
                                          ESSshifts=essS,
                                          ESSlogLik=essL))
};rm(file, shifts, shiftsM, mcmcout, postburn, essS, essL, bayes, bayesM)
dev.off()
write.table(outTable, files$outTable, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


#
#----  
#---- BAMM set final control file for all subclades, diversity estimates and proper Nshifts --------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirTrees="clades/",
             divFile="fractions/RAgamma1-1_supergroups_fractions.tsv",
             shifts=c(1, 5, 10, 50),
             dirCtlOut="bamm/controls")

if(!dir.exists(sub("\\/.*", "", files$dirCtlOut))){dir.create(sub("\\/.*", "", files$dirCtlOut))}
if(!dir.exists(files$dirCtlOut)){dir.create(files$dirCtlOut)}

div = fread(files$divFile)

trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE); i= 0
for(tree in trees){
    clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
    cat("  Working on ", clade, " (", (i = i+1), "/", length(trees), ")\n", sep="")
    
    if(!dir.exists(paste0(files$dirCtlOut, "/", clade))){dir.create(paste0(files$dirCtlOut, "/", clade))}
    
    phylo = read.tree(paste0(files$dirTrees, "/", tree))
    
    estimates = c(min(div$frac[grep(clade, div$clade)]), max(div$frac[grep(clade, div$clade)]))
    totalTips = c(min(div$total[grep(clade, div$clade)]), max(div$total[grep(clade, div$clade)]))
    
    for(j in 1:2){
        e = estimates[j]
        t = totalTips[j]
        for(s in files$shifts){
            ctl = paste0(files$dirCtlOut, "/", clade, "/control_diversification_", clade, "_div", round(e), "_shifts", s,".ctl")
            cat("   Generating", ctl, "\n")
            priors <- setBAMMpriors(phylo, total.taxa=t, outfile = NULL)
            generateControlFile(file=ctl, type = 'diversification',
                                params = list(
                                    treefile = paste0("../../../", files$dirTrees, "/", tree),
                                    outName = paste0("diversification_", clade, "_div", round(e), "_shifts", s),
                                    globalSamplingFraction = e/100,
                                    numberOfGenerations = '10000000',
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
            # system(paste("echo '\nsimulatePriorShifts = 1' >>", ctl))
        }
        
    }
}; rm(trees, i, tree, clade, phylo, estimates, totalTips, j, e, t, s, ctl, priors)


#----  
#---- BAMM analyze output interactively ------------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])

setwd("diver/Alveolata/")

files = list(tree="../../clades/clade_Alveolata.tre",
             edata="diversification_Alveolata_div31_shifts1_event_data.txt",
             mcmc= "diversification_Alveolata_div31_shifts1_mcmc_out.txt",
             shifts=1)

tree = read.tree(files$tree)

plot(tree, cex=.6)
nodelabels(frame = "none")
axisPhylo()

edata <- getEventData(tree, eventdata = files$edata, burnin=0.25)

head(edata$eventData)
bamm.tree <- plot.bammdata(edata, lwd=2, labels = T, cex = 0.5)
addBAMMshifts(bamm.tree, cex=2)
addBAMMlegend(bamm.tree)
addBAMMshifts(edata, cex=2)

edatam <- getRateThroughTimeMatrix(edata)
plotRateThroughTime(edatam, ratetype="speciation")
plotRateThroughTime(edatam, ratetype="extinction")
plotRateThroughTime(edatam, ratetype="netdiv")


# To check how many shifts are most plausible
bfmat <- computeBayesFactors(files$mcmc, expectedNumberOfShifts=files$shifts, burnin=0.1)




allrates <- getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)




mcmcout <- read.csv(files$mcmc, header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs <- table(postburn$N_shifts) / nrow(postburn)

names(post_probs)

post_probs['X'] / post_probs['Y']


shift_probs <- summary(edata)

plot.bammdata(edata, lwd=2)
plot.bammdata(edata, lwd=2, legend=T)

css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)

# plot.credibleshiftset(css)

# best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
# plot.bammdata(best, lwd = 2)
# addBAMMshifts(best, cex=2.5)

allrates <- getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="netdiv")
# plotRateThroughTime(edata, ratetype="extinction")


#----  
#---- Extract RTT table from event data file -------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])
# .rs.restartR()

# Set file names ___________________________________________________________________________________
files <- list(tree="clades/clade_Cryptista.tre",
              eventData="bamm/Cryptista/diversification_Cryptista_div44_shifts1_event_data.txt.gz")

files$clade = files$tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .) %>% sub(".*\\/", "", .)
files$shifts = files$eventData %>% sub(".*shifts", "", .) %>% sub("_.*", "", .)
files$div = files$eventData %>% sub(".*div", "", .) %>% sub("_.*", "", .)

# Open the anayzed tree ____________________________________________________________________________
tree <- read.tree(files$tree)
# plot(tree, cex=.6)
# axisPhylo()

# Read event data __________________________________________________________________________________
edata <- getEventData(tree, 
                      eventdata = files$eventData, 
                      burnin=0.25)

# head(edata$eventData)

# Plot branch rates and shifts in the tree _________________________________________________________

# http://bamm-project.org/rateshifts.html

# data(whales, events.whales)
# ed <- getEventData(whales, events.whales, burnin=0.1)
# best <- getBestShiftConfiguration(ed, expectedNumberOfShifts=1, threshold=5)
# plot.bammdata(best, lwd=1.25)
# addBAMMshifts(best, cex=2)

# library(BAMMtools)
# data(primates, events.primates)
# ed <- getEventData(primates, events.primates, burnin=0.1, type='trait')
# msc_tree <- maximumShiftCredibility(ed)

# getBranchShiftPriors

pdf(sub("_event_data.txt.gz", "_tree.pdf", files$eventData), width=11.69, height=8.27, paper='special')
bamm.tree = plot.bammdata(edata, lwd=2, labels=TRUE, cex=0.2, logcolor=TRUE)
addBAMMshifts(edata, cex=2)
addBAMMlegend(bamm.tree)
title(main = files$clade,
      sub = paste0("Speciation-extinction analysis with an estimated ", files$div, "% diversity and ", files$shifts, " expected shifts."))
dev.off()

# Now extract Rates Through time ___________________________________________________________________
edatam <- getRateThroughTimeMatrix(edata)
# or set times manually
# start.time = -ltt.plot.coords(tree)[1]
# edatam <- getRateThroughTimeMatrix(edata, start.time=start.time, end.time=0)

# quickly plot the rates
plotRateThroughTime(edatam, ratetype="netdiv")
plotRateThroughTime(edatam)
plotRateThroughTime(edatam, ratetype="extinction")

# Create a table for 'speciation' data
tmp <- as.data.frame(t(edatam$lambda))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_spe <- cbind(time=rev(edatam$times), tmp)

# Create a table for 'extinction' data
tmp <- as.data.frame(t(edatam$mu))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_ext <- cbind(time=rev(edatam$times), tmp)

# Create a table for 'diversification' data
tmp <- as.data.frame(t(edatam$lambda)) - as.data.frame(t(edatam$mu))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_div <- cbind(time=rev(edatam$times), tmp)

# Create a summary file to be exported _____________________________________________________________
# first summarise 'speciation' data
hpd = apply(edatam_spe[,-1], 1, hdi)
edatamOut = data.frame(rate="speciation",
                       time=edatam_spe$time,
                       median=apply(edatam_spe[,-1], 1, median),
                       HPDlow=hpd[1,],
                       HPDupp=hpd[2,],
                       mean=apply(edatam_spe[,-1], 1, mean),
                       sd=apply(edatam_spe[,-1], 1, sd))

# now summarise 'extinction' data
hpd = apply(edatam_ext[,-1], 1, hdi)
edatamOut = rbind(edatamOut, data.frame(rate="extinction",
                                        time=edatam_ext$time,
                                        median=apply(edatam_ext[,-1], 1, median),
                                        HPDlow=hpd[1,],
                                        HPDupp=hpd[2,],
                                        mean=apply(edatam_ext[,-1], 1, mean),
                                        sd=apply(edatam_ext[,-1], 1, sd)))

# and finally summarise 'diversification' data
hpd = apply(edatam_div[,-1], 1, hdi)
edatamOut = rbind(edatamOut, data.frame(rate="diversification",
                                        time=edatam_div$time,
                                        median=apply(edatam_div[,-1], 1, median),
                                        HPDlow=hpd[1,],
                                        HPDupp=hpd[2,],
                                        mean=apply(edatam_div[,-1], 1, mean),
                                        sd=apply(edatam_div[,-1], 1, sd)))
rm(hpd, tmp)

# Exporting the table ______________________________________________________________________________
write.table(edatamOut, sub("_event_data.txt.gz", ".tsv", files$eventData), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

#

#---- RTT BAMM output of single file ---------------------------------------------------------------

# Set file names ___________________________________________________________________________________
files <- list(rttTable="bamm/Stramenopila/diversification_Stramenopila_div34_shifts10_RTT.tsv")

files$clade = files$rttTable %>% sub("\\/diversification.*", "", .) %>% sub(".*\\/", "", .)
files$shifts = files$rttTable %>% sub(".*shifts", "", .) %>% sub("_.*", "", .)
files$div = files$rttTable %>% sub(".*div", "", .) %>% sub("_.*", "", .)

data = fread(files$rttTable)

# Plotting _________________________________________________________________________________________
# Beautifying the table for plotting
data$rate = factor(data$rate, levels=c("extinction", "speciation", "diversification"))

# Plotting
(rttPlot = ggplot(data, aes(x=-time, y=median, colour=rate))+
        geom_line()+
        geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=rate, color=NULL), alpha=0.1)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                           minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
        labs(title=paste0(files$clade),
             subtitle=paste0("Speciation-extinction analysis with an estimated ", files$div, "% diversity and ", files$shifts, " expected shifts."),
             y="Rates Through Time", x="Time before present")+
        theme_bw())

# Export plot
pdf(sub("\\.tsv", ".pdf", files$rttTable), width=11.69, height=8.27, paper='special')
plot(rttPlot)
dev.off()

#

#---- RTT BAMM output of several files from the same clade -----------------------------------------

rm(list=ls()[!ls() %in% c("geo")])
# Set file names ___________________________________________________________________________________
files <- list(dirTables="bamm/diver/Metamonada/",
              mergedOutput="bamm/diver/Metamonada_RTTs.pdf")

for(file in (grep("\\.tsv",dir(files$dirTables), value=TRUE))){
    clade = file %>% sub("diversification_", "", .) %>% sub("_.*", "", .)
    shifts = file %>% sub(".*shifts", "", .) %>% sub("_.*", "", .)
    div = file %>% sub(".*div", "", .) %>% sub("_.*", "", .)
    # Reading the table
    data = fread(paste0(files$dirTables, file))
    # Beautifying the table for plotting
    data$rate = factor(data$rate, levels=c("extinction", "speciation", "diversification"))
    # Plotting
    (rttPlot = ggplot(data, aes(x=-time, y=median, colour=rate))+
            geom_line()+
            geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=rate, color=NULL), alpha=0.1)+
            # scale_y_log10() + annotation_logticks(sides = 'l')+
            scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                               minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
            labs(title=paste0(clade),
                 subtitle=paste0("Speciation-extinction analysis with an estimated ", div, "% diversity and ", shifts, " expected shifts."),
                 y="Rates Through Time", x="Time before present")+
            theme_bw())
    # Export plot
    pdf(paste0(files$dirTables, sub("\\.tsv", ".pdf", file)), width=11.69, height=8.27, paper='special')
    plot(rttPlot)
    dev.off()
}; rm(file, data, clade, div, shifts, rttPlot)

# Pool all pdfs together
system(paste0("pdftk ", paste0(files$dirTables, grep("_RTT\\.pdf",dir(files$dirTables), value=TRUE), collapse=" "), " cat output ", files$mergedOutput))

#

#----  
#---- Plot different RTT BAMM slopes ---------------------------------------------------------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

files = list(dirFiles="diver",
             div="Min",
             shifts="50")

(files$outPlot = paste0("diver/all_RTT_div", files$div,"Estimate_shifts", files$shifts,".pdf"))

# Select all files
fileList = grep("RTT\\.tsv", dir(files$dirFiles, recursive=TRUE), value=TRUE)

# Select only files matching the number of shifts
fileList = grep(paste0("_shifts", files$shifts,"_RTT"), fileList, value=TRUE)

# And now select only those with the chosen diversity estimate
clades = unique(sub("/.*", "", fileList))
files$files = c()
for(clade in clades){
    file = grep(clade, fileList, value=TRUE)
    div = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
    if(files$div=='min' | files$div=='Min'){
        files$files = c(files$files, paste0(files$dirFiles, "/", grep(paste0("_div", min(div)), file, value=TRUE)))
    }else if(files$div=='max' | files$div=='Max'){
        files$files = c(files$files, paste0(files$dirFiles, "/", grep(paste0("_div", max(div)), file, value=TRUE)))
    }else{
        cat("\nPlease select 'files$div' as either 'min' or 'max'\n\n")
    }
}; rm(clades, clade, file, fileList, div)


data <- data.frame()
for(file in files$files){
    clade = file %>% sub("diver\\/", "", .) %>% sub("\\/.*", "", .)
    tmp = fread(file)
    tmp$clade = clade
    data = rbind(data, tmp)
};rm(file, tmp)


data$rate = factor(data$rate, levels=c("extinction", "speciation", "diversification"))

{data$colour = data$clade
    data$colour[which(data$colour=="All")]="black"
    data$colour[which(data$colour=="Amoebozoa")]="royalblue1"
    # data$colour[which(data$colour=="Breviatea")]="lightskyblue2"
    data$colour[which(data$colour=="Nucletmycea")]="steelblue2"
    data$colour[which(data$colour=="Holozoa")]="steelblue4"
    # data$colour[which(data$colour=="Ancyromonadida")]="grey20"
    # data$colour[which(data$colour=="CRuMs")]="aquamarine1"
    data$colour[which(data$colour=="Metamonada")]="forestgreen"
    data$colour[which(data$colour=="Discoba")]="orange2"
    data$colour[which(data$colour=="Haptista")]="yellow1"
    data$colour[which(data$colour=="Cryptista")]="hotpink1"
    # data$colour[which(data$colour=="Hemimastigophora")]="grey80"
    data$colour[which(data$colour=="Archaeplastida")]="darkseagreen3"
    # data$colour[which(data$colour=="Telonemia")]="plum4"
    data$colour[which(data$colour=="Rhizaria")]="darkorchid2"
    data$colour[which(data$colour=="Stramenopila")]="darkorchid3"
    data$colour[which(data$colour=="Alveolata")]="darkorchid4"}

data$clade <- factor(data$clade, 
                     levels=c("All",
                              "Amoebozoa", "Nucletmycea", "Holozoa", 
                              "Metamonada", "Discoba", 
                              "Haptista", "Cryptista", "Archaeplastida", 
                              "Rhizaria", "Stramenopila", "Alveolata"))

data$colour <- factor(data$colour, 
                      levels=c("black",
                               "royalblue1", "lightskyblue1", "steelblue2", "steelblue4",
                               "forestgreen", "orange2",
                               "yellow1", "hotpink1", "darkseagreen3",
                               "darkorchid2", "darkorchid3", "darkorchid4"))

# (outPlot = ggplot(data, aes(x=-time, y=median, colour=clade))+
#         geom_line()+
#         geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=clade, color=NULL), alpha=0.1)+
#         facet_wrap(~rate, nrow=3, scales="free")+
#         scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
#                            minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
#         labs(title=paste0("Speciation-extinction analysis with a ", files$div, " estimated diversity and ", files$shifts, " expected shifts."),
#              y="Rates Through Time", x="Time before present")+
#         scale_color_manual(values=as.character(sort(unique(data$colour))))+
#         scale_fill_manual(values=as.character(sort(unique(data$colour))))+
#         theme_bw())
# 
# pdf(files$outPlot, width=11.69, height=8.27, paper='special')
# plot(outPlot)
# dev.off()

# Now go fancy and add your subsets

datas = subset(data, clade=="Rhizaria" | clade=="Holozoa" | clade=="Archaeplastida")

(outPlots = ggplot(datas, aes(x=-time, y=median, colour=clade))+
        geom_line()+
        geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=clade, color=NULL), alpha=0.1)+
        facet_wrap(~rate, nrow=3, scales="free")+
        scale_x_continuous(breaks=seq(-max(round(datas$time, -2)), 0, 100),
                           minor_breaks=seq(-max(round(datas$time, -2)), 0, 50)) +
        labs(title=paste0("Speciation-extinction analysis with a ", files$div, " estimated diversity and ", files$shifts, " expected shifts."),
             y="Rates Through Time", x="Time before present")+
        scale_color_manual(values=as.character(sort(unique(datas$colour))))+
        scale_fill_manual(values=as.character(sort(unique(datas$colour))))+
        theme_bw())

pdf(sub("\\.pdf", "_Holo-Rhiz-Arch.pdf", files$outPlot), width=11.69, height=8.27, paper='special')
plot(outPlots)
dev.off()




# And once all PDFs are exported

system(paste0("pdftk ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\.pdf", dir(files$dirFiles), value=TRUE), collapse=" "),
              " cat output ", files$dirFiles, "/", "allClades_RTT_slopes.pdf", collapse=""))
system(paste("rm -f ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\.pdf", dir(files$dirFiles), value=TRUE), collapse=" ")))

system(paste0("pdftk ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\_Holo-Rhiz-Arch.pdf", dir(files$dirFiles), value=TRUE), collapse=" "),
              " cat output ", files$dirFiles, "/", "allClades_RTT_slopes_HoloRhizArch.pdf", collapse=""))
system(paste("rm -f ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\_Holo-Rhiz-Arch.pdf", dir(files$dirFiles), value=TRUE), collapse=" ")))


#

#---- Plot selected RTT BAMM slopes ----------------------------------------------------------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

files = list(files=c("diver/Alveolata/diversification_Alveolata_div",
    "diver/Amoebozoa/diversification_Amoebozoa_div",
    "diver/Archaeplastida/diversification_Archaeplastida_div70_shifts1_RTT.pdf",
    "diver/Cryptista/diversification_Cryptista_div",
    "diver/Discoba/diversification_Discoba_div",
    "diver/Haptista/diversification_Haptista_div",
    "diver/Holozoa/diversification_Holozoa_div",
    "diver/Metamonada/diversification_Metamonada_div",
    "diver/Nucletmycea/diversification_Nucletmycea_div",
    "diver/Rhizaria/diversification_Rhizaria_div",
    "diver/Stramenopila/diversification_Stramenopila_div"
    ),
             div="max",
             shifts="5",
             outPlot="diver/RTT_divMax_shifts5_Holo-Rhiz.pdf")

# Read  files
data <- data.frame()
for(file in files$files){
    clade = file %>% sub("diver\\/", "", .) %>% sub("\\/.*", "", .)
    tmp = fread(file)
    tmp$clade = clade
    data = rbind(data, tmp)
};rm(file, tmp)


data$rate = factor(data$rate, levels=c("extinction", "speciation", "diversification"))

{data$colour = data$clade
    data$colour[which(data$colour=="Amoebozoa")]="royalblue1"
    data$colour[which(data$colour=="Nucletmycea")]="steelblue2"
    data$colour[which(data$colour=="Holozoa")]="steelblue4"
    data$colour[which(data$colour=="Metamonada")]="forestgreen"
    data$colour[which(data$colour=="Discoba")]="orange2"
    data$colour[which(data$colour=="Haptista")]="yellow1"
    data$colour[which(data$colour=="Cryptista")]="hotpink1"
    data$colour[which(data$colour=="Archaeplastida")]="darkseagreen3"
    data$colour[which(data$colour=="Rhizaria")]="darkorchid2"
    data$colour[which(data$colour=="Stramenopila")]="darkorchid3"
    data$colour[which(data$colour=="Alveolata")]="darkorchid4"}

data$clade <- factor(data$clade, 
                     levels=c(# "Amoebozoa", 
                              # "Nucletmycea", 
                              "Holozoa", 
                              # "Metamonada", 
                              # "Discoba", 
                              # "Haptista", 
                              # "Cryptista", 
                              # "Archaeplastida", 
                              "Rhizaria"
                              # "Stramenopila", 
                              # "Alveolata"
                              ))

data$colour <- factor(data$colour, 
                      levels=c(# "royalblue1",    # Amoebozoa
                               # "steelblue2",    # Nucletmycea
                               "steelblue4",    # Holozoa
                               # "forestgreen",   # Metamonada
                               # "orange2",       # Discoba
                               # "yellow1",       # Haptista
                               # "hotpink1",      # Cryptista
                               # "darkseagreen3", # Archaeplastida
                               "darkorchid2"   # Rhizaria
                               # "darkorchid3",   # Stramenopila
                               # "darkorchid4"    # Alveolata
                               ))

(outPlot = ggplot(data, aes(x=-time, y=median, colour=clade))+
        geom_line()+
        geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=clade, color=NULL), alpha=0.1)+
        facet_wrap(~rate, nrow=3, scales="free")+
        scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                           minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
        labs(title=paste0("Speciation-extinction analysis with a ", files$div, " estimated diversity and ", files$shifts, " expected shifts."),
             y="Rates Through Time", x="Time before present")+
        scale_color_manual(values=as.character(sort(unique(data$colour))))+
        scale_fill_manual(values=as.character(sort(unique(data$colour))))+
        theme_bw())

pdf(files$outPlot, width=11.69, height=8.27, paper='special')
plot(outPlot)
dev.off()




# And once all PDFs are exported

system(paste0("pdftk ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\.pdf", dir(files$dirFiles), value=TRUE), collapse=" "),
              " cat output ", files$dirFiles, "/", "all_RTT_minDivEst.pdf", collapse=""))
system(paste("rm -f ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\.pdf", dir(files$dirFiles), value=TRUE), collapse=" ")))

system(paste0("pdftk ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\_HoloVsRhiz.pdf", dir(files$dirFiles), value=TRUE), collapse=" "),
              " cat output ", files$dirFiles, "/", "all_RTT_minDivEst_HoloVsRhiz.pdf", collapse=""))
system(paste("rm -f ", paste0(files$dirFiles, "/", grep("all_RTT.*shifts\\d+\\_HoloVsRhiz.pdf", dir(files$dirFiles), value=TRUE), collapse=" ")))


#

#---- Plot different RTT BAMM slopes OLD -----------------------------------------------------------

# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

# files <- list(files=c("tmp/Holozoa_diversity50_rate_speciation.tsv",
#                       "tmp/Metazoa_diversity50_rate_speciation.tsv",
#                       "Retaria_diversity50_rate_speciation.tsv",
#                       "Rhizaria_diversity50_rate_speciation.tsv"),
#               groups=c("Holozoa",
#                        "Metazoa",
#                        "Retaria",
#                        "Rhizaria"),
#               rate="speciation")

files <- list(files=c("",
                      "",
                      ""),
              groups=c("Div", "spe", "ext"),
              rate="Reataria")

if(length(files$files)!=length(files$group)){cat("Warning! You gave", length(files$files), "files and", length(files$groups), "group names\nPlease, give the same number of groups as input files.")}

data <- data.frame()
for(i in 1:length(files$files)){
    cat("\rWorking on '", files$groups[i], "' (", i, "/", length(files$groups), ")                    ", sep="", end="")
    tmp <- fread(files$files[i])
    median <- apply(tmp[,-1], 1, function(x) quantile(x, 0.5))
    mean <- apply(tmp[,-1], 1, mean)
    hpd <- as.data.frame(t(apply(tmp[,-1], 1, hdi)))
    tmp1 <- data.frame(group=rep(files$groups[i], nrow(tmp)),
                       time=tmp[,1],
                       median=median,
                       mean=mean,
                       hpd95=hpd$upper,
                       hpd05=hpd$lower)
    data <- rbind(data, tmp1)
};rm(i, tmp, median, mean, hpd, tmp1); cat("\rDone                                        ", end="")

# Beautifying the table ____________________________________________________________________________

data$time <- data$time * -1

geos <- subset(geo, time > min(data$time))

data$group <- factor(data$group, levels=files$groups)

# Holozoa, Metazoa, Rhizaria, Retaria
# unique(data$group)
# c("steelblue4", "steelblue2", "darkorchid2", "plum4")


# And plotting _____________________________________________________________________________________

(rttplot <- ggplot(data, aes(x=time, color=group))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_ribbon(aes(ymin=hpd05, ymax=hpd95, fill=group, color=NULL), alpha=0.1)+
        geom_line(aes(y=median), size=1) +
        # geom_line(aes(y=hpd95), alpha=0.2) +
        # geom_line(aes(y=hpd05), alpha=0.2) +
        scale_x_continuous(breaks=seq(min(round(data$time, -2)), 0, 100),
                           minor_breaks=seq(min(round(data$time, -2)), 0, 50)) +
        # scale_color_manual(values=c("steelblue4", "steelblue2", "darkorchid2", "plum4"))+
        # scale_fill_manual(values=c("steelblue4", "steelblue2", "darkorchid2", "plum4"))+
        scale_color_manual(values=c("darkorchid4", "darkorchid2", "plum4"))+
        scale_fill_manual(values=c("darkorchid4", "darkorchid2", "plum4"))+
        theme_classic()+
        labs(y=files$rate, x="Time (Ma)"))

pdf(paste0("RTT_", files$rate, "_", paste(files$groups, collapse="-"), ".pdf"), width=11.69, height=8.27, paper='special')
plot(rttplot)
dev.off()



#----  
#----  
#---- Plot the slope with the increment change on the RTT slope ------------------------------------

files = list(file="clads/Cryptista/clade_Cryptista_div44_clads_RTT.tsv",
             outPlot="clads/changes_RTT_TMP.pdf")

# Read all files, and estimate the change in diversification rate
data = fread(files$file)
tmp1 = c(); tmp2 = c()
for(i in 2:nrow(data)){
    tmp2 = c(tmp2, (data$rate[i]-data$rate[(i-1)]))
}; rm(i)
data$change = c(0, tmp2); rm(tmp2)
data$timeRev = -rev(data$time)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))

# Add a colour whether the change is positive or negative
data$colour = ifelse(data$changeSQRT < 0, "negative", "positive")
data$colour = factor(data$colour, levels=c("negative", "positive"))
data$colourSlope = c(data$colour[-1], data$colour[1])
data$colourSlope = factor(data$colourSlope, levels=c("negative", "positive"))

# And now plot
geos <- subset(geo, time>-max(data$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(barPlot = ggplot(data)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                  size=2, lineend = "round")+
        geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                 stat="identity", position=position_dodge(), alpha=0.8)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
        theme_classic()+
        scale_color_manual(values=c("orangered4", "springgreen4")) +
        scale_fill_manual(values=c("orangered4", "springgreen4")) +
        theme(legend.position="none")+
        labs(title=files$file,
             y="Change in diversification rate (square-rooted)", x="Time (Ma)"))

pdf(files$outPlot, width=11.69, height=8.27, paper='special')
plot(barPlot)
dev.off()
#

#---- Plot the slope with the increment change on the RTT slope of different files -----------------

files = list(filesDir="clads/",
             outPlot="clads/changes_RTT.pdf")
files$files = grep("RTT\\.tsv", dir(files$filesDir, recursive=TRUE), value=TRUE)

# Read all files, estimate the change in diversification rate and plot
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
for(file in files$files){
    cat("  Working on '", file, "' (", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    data = fread(paste0(files$filesDir, "/", file))
    tmp1 = c()
    for(i in 2:nrow(data)){
        tmp1 = c(tmp1, (data$rate[i]-data$rate[(i-1)]))
    }
    data$change = c(0, tmp1)
    data$timeRev = -rev(data$time)
    data$file = file
    
    # Square root the change of rates and keep whether it is positive or negative
    data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))
    
    # Add a colour whether the change is positive or negative
    data$colour = ifelse(data$changeSQRT > 0, "positive", ifelse(data$changeSQRT == 0, "zero", "negative"))
    data$colour = factor(data$colour, levels=c("positive", "zero", "negative"))
    data$colourSlope = c(data$colour[-1], data$colour[1])
    data$colourSlope = factor(data$colourSlope, levels=c("positive", "zero", "negative"))
    
    # And now plot
    geos <- subset(geo, time>-max(data$time))
    geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
    
    (barPlot = ggplot(data)+
            geom_vline(xintercept = geos$time, color="lightgrey") +
            geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                      size=2, lineend = "round")+
            geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                     stat="identity", position=position_dodge(), alpha=0.8)+
            # scale_y_log10() + annotation_logticks(sides = 'l')+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
            theme_classic()+
            scale_color_manual(values=c("springgreen4", "black", "orangered4")) +
            scale_fill_manual(values=c("springgreen4", "black", "orangered4")) +
            theme(legend.position="none")+
            labs(title=file,
                 y="Increment in diversification rate (square-rooted)", x="Time (Ma)"))
    plot(barPlot)
}; rm(file, data, tmp1, i); dev.off()

#

#----
#---- Analyze when the diversification shifts tend to happen regarding the root of the tree --------

#----  
##### Analyse ClaDS output #########################################################################
#----
#---- Export diversities for clads jobs ------------------------------------------------------------

fractions = fread("fractions/RAgamma1-1_supergroups_fractions.tsv")

for(i in unique(fractions$clade)){
    ss = subset(fractions, clade==i)
    cat(paste0('TREE="clade_', i, '.tre"\n'))
    cat(paste0('FRACTIONS="', max(ss$frac)/100, " ", min(ss$frac)/100, '"\n'))
}

tmp = data.frame()
for(i in unique(fractions$clade)){
    ss = subset(fractions, clade==i)
    cat(paste0('clade_', i, '.tre\t', max(ss$frac)/100, '\n'))
    cat(paste0('clade_', i, '.tre\t', min(ss$frac)/100, '\n'))
    tmp = rbind(tmp, data.frame(tree=paste0('clade_', i, '.tre'), fraction=max(ss$frac)/100))
    tmp = rbind(tmp, data.frame(tree=paste0('clade_', i, '.tre'), fraction=min(ss$frac)/100))
}; rm(ss, i)
write.table(tmp, "clads/fractions.tsv", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)


#----  
#---- Read Rdata from ClaDS function (julia) -------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])
# .rs.restartR()

files = list(file="clads/Nucletmycea/clade_Nucletmycea_div84_clads.Rdata",
             burn=25) # in percentage

load(files$file)


# Checking output from ClaDS Julia _________________________________________________________________

# plot(tree)
# axisPhylo()

# Plot rates of diversification in the tree
tree = CladsOutput$tree
rates = CladsOutput$lambdai_map

pdf(paste0(gsub("\\.[^\\.]+$", "_RTTtree.pdf", files$file)), width=11.69, height=8.27, paper='special')
plot_ClaDS_phylo(tree, rates, show.tip.label=TRUE, cex=0.2)
title(paste0(files$file %>% gsub("clade_", "", .) %>% gsub("_clads.*", "", .)))
dev.off()

# ltt <- as.data.frame(ltt.plot.coords(CladsOutput$tree))
# ltt$logN <- log(ltt$N, exp(1))

# Plotting the rate through time plot ______________________________________________________________

time <- apply(data.frame(CladsOutput$time_points[-length(CladsOutput$time_points)], CladsOutput$time_points[-1]),1 , mean)
clads_rtt <- data.frame(time=time, 
                        rate=CladsOutput$RTT_map)

chains <- data.frame(time=time)
burning = length(CladsOutput$rtt_chains[[1]])*(files$burn/100)
c <- 0
for(chain in CladsOutput$rtt_chains){
    c <- c + 1
    chaini <- as.data.frame(chain)
    colnames(chaini) <- paste0("Chain", c, "_iter", 1:ncol(chaini))
    chaini <- chaini[,-c(1:burning)]
    chains <- cbind(chains, chaini)
}; rm(c, chaini)
chains$toRemove <- NULL

library(HDInterval)
hpd <- as.data.frame(t(apply(t(chains), 2, hdi)))
clads_rtt$HPD_05 <- hpd$lower
clads_rtt$HPD_95 <- hpd$upper
clads_rtt$mean = apply(t(chains), 2, mean)

geos <- subset(geo, time>-max(clads_rtt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(rttPlot <-ggplot(clads_rtt, aes(x=-rev(time), y=rate))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(size=2)+
        geom_line(aes(y=HPD_05), colour="lightblue")+
        geom_line(aes(y=HPD_95), colour="lightblue")+
        # geom_line(aes(y=mean), colour="springgreen2")+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        scale_x_continuous(breaks=seq((round(max(clads_rtt$time)*-100, -2)), 0, 100),
                           minor_breaks=seq((round(min(clads_rtt$time)*-100, -2)), 0, 50)) +
        theme_classic()+
        labs(title=paste0(files$file %>% gsub("clade_", "", .) %>% gsub("_clads.*", "", .)),
             y="Diversification rate (Ln)", x="Time (Ma)"))

pdf(gsub("\\.[^\\.]+$", "_RTT.pdf", files$file), width=11.69, height=8.27, paper='special')
plot(rttPlot)
dev.off()
write.table(clads_rtt, gsub("\\.Rdata", "_RTT.tsv", files$file), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Plotting all chains for visualization purposes ___________________________________________________

# chain = melt(as.data.table(chains), id.vars="time")
# ggplot(chain, aes(x=-rev(time), y=value, colour=variable))+
#     geom_vline(xintercept = geos$time, color="lightgrey") +
#     geom_line(size=0.2)+
#     scale_y_log10() + annotation_logticks(sides = 'l')+
#     scale_x_continuous(breaks=seq((round(max(chain$time)*-100, -2)), 0, 100),
#                        minor_breaks=seq((round(min(chain$time)*-100, -2)), 0, 50)) +
#     theme_classic()+
#     theme(legend.position="none")+
#     labs(title=paste0(files$file %>% gsub("clade_", "", .) %>% gsub("_clads.*", "", .)),
#          y="Diversification rate (Ln)", x="Time (Ma)")


# Plotting the Lineages through time plot __________________________________________________________

clads_ltt <- data.frame(time_points=CladsOutput$time_points,
                        time=-rev(CladsOutput$time_points),
                        lineages=CladsOutput$DTT_mean)

enhanced = data.frame(time=-rev(CladsOutput$time_points))
for(t in 1:length(CladsOutput$enhanced_trees)){
    cat("\r  Working on enhanced tree ", t, "/", length(CladsOutput$enhanced_trees), sep="")
    tmp = as.data.frame(ltt.plot.coords(CladsOutput$enhanced_trees[[t]]$tree))
    tmp1 = c()
    for(ti in enhanced$time){
        i = ifelse(ti < min(tmp$time), 1, max(subset(tmp, time < ti)$N)+1)
        tmp1 = c(tmp1, i)
    }
    enhanced = cbind(enhanced, tmp1)
}; rm(t, tmp, tmp1, ti); cat("\n")
colnames(enhanced) = c("time", paste0("tree", 1:(ncol(enhanced)-1)))

clads_ltt$mean = apply(enhanced[,-1], 1, mean)
clads_ltt$min = apply(enhanced[,-1], 1, min)
clads_ltt$max = apply(enhanced[,-1], 1, max)

geos <- subset(geo, time>min(clads_ltt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(lttPlot <-ggplot(clads_ltt, aes(x=time, y=lineages))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        # geom_point()+
        # geom_line(aes(y=min), colour="lightblue")+
        # geom_line(aes(y=max), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree1), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree2), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree3), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree4), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree5), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree6), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree7), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree8), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree9), colour="lightblue")+
        geom_line(data=enhanced, aes(x=time, y=tree10), colour="lightblue")+
        geom_line(size=2)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        scale_x_continuous(breaks=seq((round(min(clads_ltt$time), -2)), 0, 100),
                           minor_breaks=seq((round(min(clads_ltt$time), -2)), 0, 50)) +
        theme_classic()+
        labs(title=paste0(files$file %>% gsub("clade_", "", .) %>% gsub("_clads.*", "", .)),
             y="Diversity Through Time (Ln)", x="Time (Ma)"))

pdf(gsub("\\.[^\\.]+$", "_DTT.pdf", files$file), width=11.69, height=8.27, paper='special')
plot(lttPlot)
dev.off()
write.table(clads_ltt, gsub("\\.Rdata", "_DTT.tsv", files$file), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


#

#---- Read Rdata from several ClaDS files (julia) --------------------------------------------------

# Checking output from two different ClaDS Julia output ____________________________________________
# .rs.restartR()
setwd("clads/")
rm(list=ls()[!ls() %in% c("geo")])

groups <- list(Amoebozoa="Amoebozoa/clade_Amoebozoa_div89_clads.Rdata",
               #Nucletmycea,
               Holozoa="Holozoa/clade_Holozoa_div83_clads.Rdata",
               Discoba="Discoba/clade_Discoba_div84_clads.Rdata",
               Metamonada="Metamonada/clade_Metamonada_div81_clads.Rdata",
               Haptista="Haptista/clade_Haptista_div80_clads.Rdata",
               Cryptista="Cryptista/clade_Cryptista_div44_clads.Rdata",
               Archaeplastida="Archaeplastida/clade_Archaeplastida_div70_clads.Rdata",
               Rhizaria="Rhizaria/clade_Rhizaria_div72_clads.Rdata",
               Stramenopila="Stramenopila/clade_Stramenopila_div70_clads.Rdata"
               #Alveolata
)

getRTTfromCladsOutput <- function(cladsOutput, verbose=TRUE){
    if(verbose){cat("  Reading file\n")}
    tmp <- apply(data.frame(CladsOutput$time_points[-length(CladsOutput$time_points)], CladsOutput$time_points[-1]),1 , mean)
    clads_rtt <- data.frame(time=tmp, 
                            rate=CladsOutput$RTT_map); rm(tmp)
    
    if(verbose){cat("  Extracting chains\n")}
    burn <- 100
    chains <- data.frame(toRemove=rep(NA, length(cladsOutput$rtt_chains[[1]][[1]])))
    c <- 0
    for(chain in cladsOutput$rtt_chains){
        c <- c + 1
        cat("\r    ", c, "/", length(cladsOutput$rtt_chains), sep="")
        chaini <- as.data.frame(chain)
        colnames(chaini) <- paste0("Chain", c, "_iter", 1:ncol(chaini))
        chaini <- chaini[,-c(1:burn)]
        chains <- cbind(chains, chaini)
    }; rm(c, chaini)
    chains$toRemove <- NULL
    
    if(verbose){cat("\n  Calculating HPD intervals\n")}
    hpd <- as.data.frame(t(apply(t(chains), 2, hdi)))
    clads_rtt$HPD_05 <- hpd$lower
    clads_rtt$HPD_95 <- hpd$upper
    
    clads_rtt$revTime <- -rev(clads_rtt$time)
    if(verbose){cat("Done\n")}
    return(clads_rtt)
}

rtt = data.frame()
dtt = data.frame()
for(g in names(groups)){
    cat("Loading data for '", g ,"' (", grep(g, names(groups)), "/", length(groups), ")\n", sep="")
    load(groups[[g]])
    tmp <- CladsOutput
    
    tmpt <- getRTTfromCladsOutput(cladsOutput=tmp)
    tmpt$group <- g
    
    dtt = rbind(dtt, data.frame(time_points=CladsOutput$time_points,
                                time=-rev(CladsOutput$time_points),
                                lineages=CladsOutput$DTT_mean,
                                group=g))
    
    rtt <- rbind(rtt, tmpt)
}; rm(g, CladsOutput, tmp, tmpt)


geos <- subset(geo, time>-max(rtt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

{rtt$colour = rtt$group
    rtt$colour[which(rtt$colour=="Main_tree")]="black"
    rtt$colour[which(rtt$colour=="Amoebozoa")]="royalblue1"
    rtt$colour[which(rtt$colour=="Nucletmycea")]="steelblue2"
    rtt$colour[which(rtt$colour=="Holozoa")]="steelblue4"
    rtt$colour[which(rtt$colour=="Metamonada")]="forestgreen"
    rtt$colour[which(rtt$colour=="Discoba")]="orange2"
    rtt$colour[which(rtt$colour=="Haptista")]="yellow1"
    rtt$colour[which(rtt$colour=="Cryptista")]="hotpink1"
    rtt$colour[which(rtt$colour=="Archaeplastida")]="darkseagreen3"
    rtt$colour[which(rtt$colour=="Rhizaria")]="darkorchid2"
    rtt$colour[which(rtt$colour=="Stramenopila")]="darkorchid3"
    rtt$colour[which(rtt$colour=="Alveolata")]="darkorchid4"}

rtt$group = factor(rtt$group, levels=c(#"All",
    "Amoebozoa", #"Nucletmycea", "Holozoa",
    "Metamonada", "Discoba",
    "Haptista", "Cryptista", "Archaeplastida",
    "Rhizaria", "Stramenopila", "Alveolata"))

rtt$colour <- factor(rtt$colour,
                     levels=c(#"black",
                         "royalblue1", #"steelblue2", "steelblue4",
                         "forestgreen", "orange2",
                         "yellow1", "hotpink1", "darkseagreen3",
                         "darkorchid2", "darkorchid3", "darkorchid4"))

(rttplot <- ggplot(rtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=revTime, y=rate, colour=group), size=1)+
        geom_ribbon(aes(x=revTime, ymax=HPD_95, ymin=HPD_05, fill=group), alpha=0.2, colour = NA) +
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.01, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.005, label=period), size=2)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(rtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(rtt$colour))))+
        labs(y="Diversification rate (Ln)", x="Time (Ma)"))

pdf("clads_RTT_all.pdf", width=11.69, height=8.27, paper='special')
plot(rttplot)
dev.off()


{dtt$colour = dtt$group
    dtt$colour[which(dtt$colour=="Main_tree")]="black"
    dtt$colour[which(dtt$colour=="Amoebozoa")]="royalblue1"
    dtt$colour[which(dtt$colour=="Nucletmycea")]="steelblue2"
    dtt$colour[which(dtt$colour=="Holozoa")]="steelblue4"
    dtt$colour[which(dtt$colour=="Metamonada")]="forestgreen"
    dtt$colour[which(dtt$colour=="Discoba")]="orange2"
    dtt$colour[which(dtt$colour=="Haptista")]="yellow1"
    dtt$colour[which(dtt$colour=="Cryptista")]="hotpink1"
    dtt$colour[which(dtt$colour=="Archaeplastida")]="darkseagreen3"
    dtt$colour[which(dtt$colour=="Rhizaria")]="darkorchid2"
    dtt$colour[which(dtt$colour=="Stramenopila")]="darkorchid3"
    dtt$colour[which(dtt$colour=="Alveolata")]="darkorchid4"}

dtt$group = factor(dtt$group, levels=c(#"All",
    "Amoebozoa", #"Nucletmycea", "Holozoa",
    "Metamonada", "Discoba",
    "Haptista", "Cryptista", "Archaeplastida",
    "Rhizaria", "Stramenopila", "Alveolata"))

dtt$colour <- factor(dtt$colour,
                     levels=c(#"black",
                         "royalblue1", #"steelblue2", "steelblue4",
                         "forestgreen", "orange2",
                         "yellow1", "hotpink1", "darkseagreen3",
                         "darkorchid2", "darkorchid3", "darkorchid4"))

(dttplot <- ggplot(dtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=time, y=lineages, colour=group), size=1)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(dtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(dtt$colour))))+
        labs(y="Diversification rate (Ln)", x="Time (Ma)"))

pdf("clads_DTT_all.pdf", width=11.69, height=8.27, paper='special')
plot(dttplot)
dev.off()

#
#---- Read Rdata from several processed ClaDS tables (julia) ---------------------------------------

# Checking output from two different ClaDS Julia output ____________________________________________
# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

groups <- list(# Amoebozoa="clads/Amoebozoa/clade_Amoebozoa_div89_clads_RTT.tsv",
               # Nucletmycea="clads/Nucletmycea/clade_Nucletmycea_div84_clads_RTT.tsv",
               Holozoa="clads/Holozoa/clade_Holozoa_div83_clads_RTT.tsv",
               # Discoba="clads/Discoba/clade_Discoba_div84_clads_RTT.tsv",
               # Metamonada="clads/Metamonada/clade_Metamonada_div81_clads_RTT.tsv",
               # Haptista="clads/Haptista/clade_Haptista_div80_clads_RTT.tsv",
               # Cryptista="clads/Cryptista/clade_Cryptista_div44_clads_RTT.tsv",
               Archaeplastida="clads/Archaeplastida/clade_Archaeplastida_div70_clads_RTT.tsv",
               Rhizaria="clads/Rhizaria/clade_Rhizaria_div71_clads_RTT.tsv"
               # Stramenopila="clads/Stramenopila/clade_Stramenopila_div70_clads_RTT.tsv"
               # Alveolata="clads/Alveolata"
)

rtt = data.frame()
dtt = data.frame(); DTT=TRUE
for(g in names(groups)){
    cat("Loading data for '", g ,"' (", grep(g, names(groups)), "/", length(groups), ")\n", sep="")
    
    tmp = fread(groups[[g]])
    tmp$revTime = -rev(tmp$time)
    tmp$group = g
    rtt = rbind(rtt, tmp)
    
    if(DTT){
        tmp = fread(gsub("RTT", "DTT", groups[[g]]))
        tmp$group = g
        tmp$revTime = -rev(tmp$time)
        dtt = rbind(dtt, tmp)
    }
}; rm(g, tmp, DTT)

geos <- subset(geo, time>-max(rtt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

{rtt$colour = rtt$group
    rtt$colour[which(rtt$colour=="Main_tree")]="black"
    rtt$colour[which(rtt$colour=="Amoebozoa")]="royalblue1"
    rtt$colour[which(rtt$colour=="Nucletmycea")]="steelblue2"
    rtt$colour[which(rtt$colour=="Holozoa")]="steelblue4"
    rtt$colour[which(rtt$colour=="Metamonada")]="forestgreen"
    rtt$colour[which(rtt$colour=="Discoba")]="orange2"
    rtt$colour[which(rtt$colour=="Haptista")]="yellow1"
    rtt$colour[which(rtt$colour=="Cryptista")]="hotpink1"
    rtt$colour[which(rtt$colour=="Archaeplastida")]="darkseagreen3"
    rtt$colour[which(rtt$colour=="Rhizaria")]="darkorchid2"
    rtt$colour[which(rtt$colour=="Stramenopila")]="darkorchid3"
    rtt$colour[which(rtt$colour=="Alveolata")]="darkorchid4"}

rtt$group = factor(rtt$group, levels=c(#"All",
    "Amoebozoa", 
    "Nucletmycea",
    "Holozoa",
    "Metamonada",
    "Discoba",
    "Haptista",
    "Cryptista",
    "Archaeplastida",
    "Rhizaria",
    "Stramenopila"
    # "Alveolata"
    ))

rtt$colour <- factor(rtt$colour,
                     levels=c(#"black",
                         "royalblue1",
                         "steelblue2",
                         "steelblue4",
                         "forestgreen",
                         "orange2",
                         "yellow1",
                         "hotpink1",
                         "darkseagreen3",
                         "darkorchid2", 
                         "darkorchid3"
                         # "darkorchid4"
                         ))

(rttplot <- ggplot(rtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_ribbon(aes(x=revTime, ymax=HPD_95, ymin=HPD_05, fill=group), alpha=0.2, colour = NA) +
        geom_line(aes(x=revTime, y=rate, colour=group), size=1)+
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.01, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.005, label=period), size=2)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(rtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(rtt$colour))))+
        labs(y="Diversification rate (Ln)", x="Time (Ma)"))

pdf("clads/clads_RTT_holoRhizArch.pdf", width=11.69, height=8.27, paper='special')
plot(rttplot)
dev.off()


{dtt$colour = dtt$group
    dtt$colour[which(dtt$colour=="Main_tree")]="black"
    dtt$colour[which(dtt$colour=="Amoebozoa")]="royalblue1"
    dtt$colour[which(dtt$colour=="Nucletmycea")]="steelblue2"
    dtt$colour[which(dtt$colour=="Holozoa")]="steelblue4"
    dtt$colour[which(dtt$colour=="Metamonada")]="forestgreen"
    dtt$colour[which(dtt$colour=="Discoba")]="orange2"
    dtt$colour[which(dtt$colour=="Haptista")]="yellow1"
    dtt$colour[which(dtt$colour=="Cryptista")]="hotpink1"
    dtt$colour[which(dtt$colour=="Archaeplastida")]="darkseagreen3"
    dtt$colour[which(dtt$colour=="Rhizaria")]="darkorchid2"
    dtt$colour[which(dtt$colour=="Stramenopila")]="darkorchid3"
    dtt$colour[which(dtt$colour=="Alveolata")]="darkorchid4"}

dtt$group = factor(dtt$group, levels=c(#"All",
    "Amoebozoa",
    "Nucletmycea",
    "Holozoa",
    "Metamonada",
    "Discoba",
    "Haptista",
    "Cryptista",
    "Archaeplastida",
    "Rhizaria",
    "Stramenopila"
    # "Alveolata"
))

dtt$colour <- factor(dtt$colour,
                     levels=c(#"black",
                         "royalblue1",
                         "steelblue2",
                         "steelblue4",
                         "forestgreen",
                         "orange2",
                         "yellow1",
                         "hotpink1",
                         "darkseagreen3",
                         "darkorchid2",
                         "darkorchid3"
                         # "darkorchid4"
                     ))

(dttplot <- ggplot(dtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=time, y=lineages, colour=group), size=1)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(dtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(dtt$colour))))+
        labs(y="Estimated diversity through time (Ln)", x="Time (Ma)"))

pdf("clads/clads_DTT_holoRhizArch.pdf", width=11.69, height=8.27, paper='special')
plot(dttplot)
dev.off()

#
#----
#---- Plot the slope with the increment change on the RTT slope ------------------------------------

files = list(file="clads/Cryptista/clade_Cryptista_div44_clads_RTT.tsv",
             outPlot="clads/changes_RTT_TMP.pdf")

# Read all files, and estimate the change in diversification rate
data = fread(files$file)
tmp1 = c(); tmp2 = c()
for(i in 2:nrow(data)){
    tmp2 = c(tmp2, (data$rate[i]-data$rate[(i-1)]))
}; rm(i)
data$change = c(0, tmp2); rm(tmp2)
data$timeRev = -rev(data$time)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))

# Add a colour whether the change is positive or negative
data$colour = ifelse(data$changeSQRT < 0, "negative", "positive")
data$colour = factor(data$colour, levels=c("negative", "positive"))
data$colourSlope = c(data$colour[-1], data$colour[1])
data$colourSlope = factor(data$colourSlope, levels=c("negative", "positive"))

# And now plot
geos <- subset(geo, time>-max(data$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(barPlot = ggplot(data)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                  size=2, lineend = "round")+
        geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                 stat="identity", position=position_dodge(), alpha=0.8)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
        theme_classic()+
        scale_color_manual(values=c("orangered4", "springgreen4")) +
        scale_fill_manual(values=c("orangered4", "springgreen4")) +
        theme(legend.position="none")+
        labs(title=files$file,
             y="Change in diversification rate (square-rooted)", x="Time (Ma)"))

pdf(files$outPlot, width=11.69, height=8.27, paper='special')
plot(barPlot)
dev.off()
#

#---- Plot the slope with the increment change on the RTT slope of different files -----------------

files = list(filesDir="clads/",
             outPlot="clads/changes_RTT.pdf")
files$files = grep("RTT\\.tsv", dir(files$filesDir, recursive=TRUE), value=TRUE)

# Read all files, estimate the change in diversification rate and plot
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
for(file in files$files){
    cat("  Working on '", file, "' (", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    data = fread(paste0(files$filesDir, "/", file))
    tmp1 = c()
    for(i in 2:nrow(data)){
        tmp1 = c(tmp1, (data$rate[i]-data$rate[(i-1)]))
    }
    data$change = c(0, tmp1)
    data$timeRev = -rev(data$time)
    data$file = file
    
    # Square root the change of rates and keep whether it is positive or negative
    data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))
    
    # Add a colour whether the change is positive or negative
    data$colour = ifelse(data$changeSQRT > 0, "positive", ifelse(data$changeSQRT == 0, "zero", "negative"))
    data$colour = factor(data$colour, levels=c("positive", "zero", "negative"))
    data$colourSlope = c(data$colour[-1], data$colour[1])
    data$colourSlope = factor(data$colourSlope, levels=c("positive", "zero", "negative"))
    
    # And now plot
    geos <- subset(geo, time>-max(data$time))
    geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
    
    (barPlot = ggplot(data)+
            geom_vline(xintercept = geos$time, color="lightgrey") +
            geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                      size=2, lineend = "round")+
            geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                     stat="identity", position=position_dodge(), alpha=0.8)+
            # scale_y_log10() + annotation_logticks(sides = 'l')+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
            theme_classic()+
            scale_color_manual(values=c("springgreen4", "black", "orangered4")) +
            scale_fill_manual(values=c("springgreen4", "black", "orangered4")) +
            theme(legend.position="none")+
            labs(title=file,
                 y="Increment in diversification rate (square-rooted)", x="Time (Ma)"))
    plot(barPlot)
}; rm(file, data, tmp1, i); dev.off()

#

#----
#---- Trying to find diversification-dependent patterns --------------------------------------------

files = list(filesDir="clads/",
             outPlot="clads/changes_RTT_allVsAll.pdf")

files$files = grep("RTT\\.tsv", dir(files$filesDir, recursive=TRUE), value=TRUE)

# Read all files, and estimate the change in diversification rate
data = data.frame(); tmp1 = c(); tmp2 = c()
for(file in files$files){
    tmp = fread(paste0(files$filesDir, "/", file))
    tmp1 = c(); tmp2 = c()
    for(i in 2:nrow(tmp)){
        tmp2 = c(tmp2, (tmp$rate[i]-tmp$rate[(i-1)]))
    }
    tmp$file = file
    tmp$change = c(0, tmp2)
    tmp$timeRev = -rev(tmp$time)
    data = rbind(data, tmp)
}; rm(file, tmp, tmp1, tmp2, i)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))

# And now plot all pair-wise comparisons
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
done = c()
for(file1 in files$files){
    for( file2 in files$files){
        if(file1 != file2 & ((!paste0(file1, "-", file2) %in% done) | (!paste0(file2, "-", file1) %in% done))){
            file1="Rhizaria/clade_Rhizaria_div71_clads_RTT.tsv"
            # file2="Holozoa/clade_Holozoa_div83_clads_RTT.tsv"
            file2="Stramenopila/clade_Stramenopila_div70_clads_RTT.tsv"
            cat("  Plotting", file1, "against", file2, "\n")
            done = c(done, paste0(file1, "-", file2))
            datas = subset(data, file==file1 | file==file2)
            
            geos <- subset(geo, time>-max(datas$time))
            geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
            
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="pearson")
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="kendall")
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="spearman")
            
            (barPlot = ggplot(datas, aes(x=timeMidRev, y=changeSQRT, fill=file, group=file))+
                    geom_vline(xintercept = geos$time, color="lightgrey") +
                    geom_line(aes(x=timeRev, y=sqrt(rate), color=file))+
                    geom_bar(stat="identity", position=position_dodge(), width=10, alpha=0.8)+
                    # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
                    # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
                    theme_classic()+
                    theme(legend.position="bottom")+
                    labs(y="Change in diversification rate (square-rooted)", x="Time (Ma)"))
            plot(barPlot)
        }
    }
}; rm(done, file1, file2, tmp); dev.off()

#

#----
#----
#----
#----
#----
#---- FROM HERE ON ARE OLD SCRIPTS THAT I DON'T DARE TRASHING BUT I SHOULD -------------------------
#----
# #---- Estimate phylogenetic species richness based on the rarefied trees and abundance data --------
# 
# rm(list=ls()[!ls() %in% c("geo")])
# # .rs.restartR()
# 
# files = list(abundances="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/resources/tipNames_OTUreads.tsv",
#              phylo="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/step3f/raxml-cat/rootD/step3f_RAcat1_iqtreef_GTRg_rep1_cleaned_rootD_coloured.tre",
#              dirTrees="clades",
#              dirRare="phylo",
#              outSlopes="RAcat1-1_supergroups_rarefied_fraction.pdf",
#              outData=  "RAcat1-1_supergroups_rarefied_fraction.tsv",
#              outDir="fractions")
# 
# # Create directories if doesn't exist
# if(!dir.exists(files$outDir)){dir.create(files$outDir)}
# if(!dir.exists(files$dirTrees)){dir.create(files$dirTrees)}
# if(!dir.exists(paste0(files$dirTrees, "/", files$dirRare))){dir.create(paste0(files$dirTrees, "/", files$dirRare))}
# system(paste0("cp ", files$phylo, " ", files$dirTrees, "/", files$dirRare, "/", basename(files$phylo)))
# 
# # Read the phylogenetic tree and export subclades __________________________________________________
# 
# # treei <- read.beast(files$phylo)
# # tree <- as.phylo(treei)
# # 
# # # Transform files
# # treed <- as_tibble(treei)
# # names(treed) <- gsub("!", "", names(treed))
# # 
# # # Select all annotated nodes 
# # (annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))
# # 
# # # Extract clades with more than 300 tips 
# # i = 0
# # for(annot in annotations){
# #     cat("  Extracting ", annot, " (", (i = i + 1), "/", length(annotations), ")\n", sep="")
# #     node = as.numeric(subset(treed, name==annot)$node)
# #     tmp = tree_subset(treei, node=node, levels_back=0)
# #     tmp = as.phylo(tmp)
# #     ntips = length(tmp$tip.label)
# #     if(ntips >= 300){
# #         out = paste0(files$dirTrees, "/", files$dirRare, "/clade_", annot, ".tre")
# #         cat("    Exporting tree to:", out, "with", ntips, "tips\n")
# #         write.tree(tmp, file=out)
# #     }else{
# #         cat("    Clade", annot, "has", ntips, "tips, so it is not exported\n")
# #     }
# # };rm(i, annot, node, tmp, out, ntips)
# # system(paste0("rm -f ", files$dirTrees, "/", files$dirRare, "/", basename(files$phylo)))
# 
# # Rarefy the phylogentic trees _____________________________________________________________________
# for(tree in grep("^clade_.*\\.tre$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     system(paste0("treeRarefy.py -t ", files$dirTrees, "/", files$dirRare, "/", tree, " -a ", files$abundances, " -r 0", round(Ntip(tree)*10, -3)))
# }; rm(tree, phylo, maxRange)
# 
# # Now get the slopes _______________________________________________________________________________
# 
# pdf(files$outSlopes, width=11.69, height=8.27, paper='special')
# outData = data.frame()
# for(table in grep("^clade_.*rarefied.*\\.tsv$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     clade = table %>% sub("clade_", "", .) %>% sub("_.*", "", .)
#     cat("  Working on", clade, "\n")
#     tree = read.tree(paste0(files$dirTrees, "/", files$dirRare, "/",
#                             grep(paste0(clade, "\\.tre$"), dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)))
#     data = fread(paste0(files$dirTrees, "/", files$dirRare, "/", table))
#     
#     # Apply an Asymptotic regression model
#     # Y=a−(a−b)exp(−cX)
#     # a : maximum attainable Y
#     # b : Y at x=0 
#     # c : proportional to the relative rate of Y increase while X increases
#     model <- drm(data$mean ~ data$sampleSize, fct = DRC.asymReg())
#     plateau = summary(model)$coefficients[3]
#     data$fitted = fitted(model)
#     div = plateau - (plateau - summary(model)$coefficients[1]) * exp(-summary(model)$coefficients[2] * Ntip(tree))
#     frac = div/plateau*100
#     divt = (Ntip(tree)*100)/(frac)
#     
#     outData = rbind(outData, data.frame(clade= clade, 
#                                         tips=Ntip(tree), 
#                                         fraction=frac,
#                                         total=divt))
#     
#     slope = ggplot(data)+
#         geom_point(aes(x=sampleSize, y=mean))+
#         geom_line(aes(x=sampleSize, y=fitted), color="orangered3")+
#         geom_hline(yintercept=plateau, color="orangered3")+
#         geom_segment(aes(x=Ntip(tree), y=0, 
#                          xend=Ntip(tree), yend=div), colour="springgreen3")+
#         geom_segment(aes(x=0, y=div, 
#                          xend=Ntip(tree), yend=div), colour="springgreen3")+
#         theme_minimal()+
#         labs(title=paste0(clade, ": ", Ntip(tree), " tips"),
#              subtitle=paste0("This tree represents an estimated ", round(frac, 2), "% of the global diversity (",
#                              "estimated to be ", round(divt, 1), " total tips)."),
#              y="Mean phylogenetic distance to Root", x="Sample size")
#     slope
#     plot(slope)
# }; rm(table, clade, data, model, plateau, div, frac, divt, slope)
# dev.off()
# write.table(outData, files$outData, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
# 
# #---- Estimate phylogenetic species richness based on the rarefied trees ---------------------------
# 
# # .rs.restartR()
# rm(list=ls()[!ls() %in% c("geo")]); gc()
# 
# files = list(phylo="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/step3f/raxml-cat/rootD/step3f_RAcat1_iqtreef_GTRg_rep1_cleaned_rootD_coloured.tre",
#              dirTrees="clades",
#              dirRare="rarefy",
#              outSlopes="RAcat1-1_supergroups_rarefiedNotAbun_fraction.pdf",
#              outData=  "RAcat1-1_supergroups_rarefiedNotAbun_fraction.tsv",
#              outDir="fractions")
# 
# # Create directory if doesn't exist
# if(!dir.exists(files$outDir)){dir.create(files$outDir)}else{cat("  directory '", files$outDir, "' already exists\n", sep="")}
# 
# # if(!dir.exists(files$dirTrees)){dir.create(files$dirTrees)}
# # if(!dir.exists(paste0(files$dirTrees, "/", files$dirRare))){dir.create(paste0(files$dirTrees, "/", files$dirRare))}
# # system(paste0("cp ", files$phylo, " ", files$dirTrees, "/", files$dirRare, "/", basename(files$phylo)))
# 
# # Read the phylogenetic tree and export subclades __________________________________________________
# 
# # treei <- read.beast(files$phylo)
# # tree <- as.phylo(treei)
# # 
# # # Transform files
# # treed <- as_tibble(treei)
# # names(treed) <- gsub("!", "", names(treed))
# # 
# # # Select all annotated nodes 
# # (annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))
# # 
# # # Extract clades with more than 300 tips 
# # i = 0
# # for(annot in annotations){
# #     cat("  Extracting ", annot, " (", (i = i + 1), "/", length(annotations), ")\n", sep="")
# #     node = as.numeric(subset(treed, name==annot)$node)
# #     tmp = tree_subset(treei, node=node, levels_back=0)
# #     tmp = as.phylo(tmp)
# #     ntips = length(tmp$tip.label)
# #     if(ntips >= 300){
# #         out = paste0(files$dirTrees, "/", files$dirRare, "/clade_", annot, ".tre")
# #         cat("    Exporting tree to:", out, "with", ntips, "tips\n")
# #         write.tree(tmp, file=out)
# #     }else{
# #         cat("    Clade", annot, "has", ntips, "tips, so it is not exported\n")
# #     }
# # };rm(i, annot, node, tmp, out, ntips)
# 
# # Rarefy the trees _________________________________________________________________________________
# for(tree in grep("^clade_.*\\.tre$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     cat(" ", tree, "\n")
#     file = read.tree(paste0(files$dirTrees, "/", files$dirRare, "/", tree))
#     max = round(Ntip(file)*10, -3)
#     cat(paste0("treeRarefy.py -t ", files$dirTrees, "/", files$dirRare, "/", tree, " -r 0 ", max, " -s ", max/100), "\n")
#     # system(paste0("treeRarefy.py -t ", files$dirTrees, "/", files$dirRare, "/", tree, " -r 0 ", max, " -s ", max/100))
# }; rm(tree, phylo, maxRange, file)
# 
# # Now get the slopes _______________________________________________________________________________
# 
# pdf(files$outSlopes, width=11.69, height=8.27, paper='special')
# outData = data.frame()
# for(table in grep("^clade_.*rarefied.*\\.tsv$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     clade = table %>% sub("clade_", "", .) %>% sub("_.*", "", .)
#     cat("  Working on", clade, "\n")
#     tree = read.tree(paste0(files$dirTrees, "/", files$dirRare, "/",
#                             grep(paste0(clade, "\\.tre$"), dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)))
#     
#     data = fread(paste0(files$dirTrees, "/", files$dirRare, "/", table))
#     
#     # Estimate the sample size needed to arrive at a 99% diversity
#     tips=Ntip(tree)
#     plateau = sum(node.depth.edgelength(tree)[0:tips])
#     div99= subset(data, mean==min(data$mean[data$mean > plateau*0.99]))$sampleSize
#     frac = (tips*99)/div99
#     divt = (tips*100)/(frac)
#     
#     outData = rbind(outData, data.frame(clade= clade, 
#                                         tips=tips, 
#                                         fraction=frac,
#                                         total=divt))
#     
#     (slope = ggplot(data)+
#             geom_point(aes(x=sampleSize, y=mean))+
#             # geom_line(aes(x=sampleSize, y=fitted), color="orangered3")+
#             geom_hline(yintercept=plateau, color="orangered3")+
#             geom_vline(xintercept=div99, color="orangered3")+
#             geom_segment(aes(x=Ntip(tree), y=0, 
#                              xend=Ntip(tree), yend=div99), colour="springgreen3")+
#             geom_segment(aes(x=0, y=div99, 
#                              xend=Ntip(tree), yend=div99), colour="springgreen3")+
#             theme_minimal()+
#             labs(title=paste0(clade, ": ", Ntip(tree), " tips"),
#                  subtitle=paste0("This tree represents an estimated ~", round(frac, 2), "% of the global diversity (",
#                                  "estimated to be ", round(divt, 1), " total tips)."),
#                  y="Mean phylogenetic distance to Root", x="Sample size"))
#     plot(slope)
# }; rm(table, clade, data, plateau, div99, frac, divt, slope); dev.off()
# write.table(outData, files$outData, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
# 
# #---- Estimate phylogenetic species richness based on the rarefied trees with 'Rarefy' -------------
# 
# # .rs.restartR()
# rm(list=ls()[!ls() %in% c("geo")]); gc()
# library(Rarefy)
# 
# files = list(phylo="/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/step3f/raxml-cat/rootD/step3f_RAcat1_iqtreef_GTRg_rep1_cleaned_rootD_coloured.tre",
#              dirTrees="clades",
#              dirRare="phylo",
#              outSlopes="RAcat1-1_supergroups_rarefiedNotAbun_fraction.pdf",
#              outData=  "RAcat1-1_supergroups_rarefiedNotAbun_fraction.tsv",
#              outDir="fractions")
# 
# # Create directory if doesn't exist
# if(!dir.exists(files$outDir)){dir.create(files$outDir)}else{cat("  directory '", files$outDir, "' already exists\n", sep="")}
# 
# # Rarefy the trees _________________________________________________________________________________
# for(tree in grep("^clade_.*\\.tre$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     cat(" ", tree, "\n")
#     phylo = read.tree(paste0(files$dirTrees, "/", files$dirRare, "/", tree))
#     # plot(phylo, show.tip.label = FALSE)
#     rare = rare_alpha(node.depth.edgelength(phylo)[0:Ntip(phylo)])
#     rare = rare_phylo(tree=tree)
# }; rm(tree, phylo, maxRange, file)
# 
# 
# # Now get the slopes _______________________________________________________________________________
# 
# pdf(files$outSlopes, width=11.69, height=8.27, paper='special')
# outData = data.frame()
# for(table in grep("^clade_.*rarefied.*\\.tsv$", dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)){
#     clade = table %>% sub("clade_", "", .) %>% sub("_.*", "", .)
#     cat("  Working on", clade, "\n")
#     tree = read.tree(paste0(files$dirTrees, "/", files$dirRare, "/",
#                             grep(paste0(clade, "\\.tre$"), dir(paste0(files$dirTrees, "/", files$dirRare)), value=TRUE)))
#     
#     data = fread(paste0(files$dirTrees, "/", files$dirRare, "/", table))
#     
#     # Estimate the sample size needed to arrive at a 99% diversity
#     tips=Ntip(tree)
#     plateau = sum(node.depth.edgelength(tree)[0:tips])
#     div99= subset(data, mean==min(data$mean[data$mean > plateau*0.99]))$sampleSize
#     frac = (tips*99)/div99
#     divt = (tips*100)/(frac)
#     
#     outData = rbind(outData, data.frame(clade= clade, 
#                                         tips=tips, 
#                                         fraction=frac,
#                                         total=divt))
#     
#     (slope = ggplot(data)+
#             geom_point(aes(x=sampleSize, y=mean))+
#             # geom_line(aes(x=sampleSize, y=fitted), color="orangered3")+
#             geom_hline(yintercept=plateau, color="orangered3")+
#             geom_vline(xintercept=div99, color="orangered3")+
#             geom_segment(aes(x=Ntip(tree), y=0, 
#                              xend=Ntip(tree), yend=div99), colour="springgreen3")+
#             geom_segment(aes(x=0, y=div99, 
#                              xend=Ntip(tree), yend=div99), colour="springgreen3")+
#             theme_minimal()+
#             labs(title=paste0(clade, ": ", Ntip(tree), " tips"),
#                  subtitle=paste0("This tree represents an estimated ~", round(frac, 2), "% of the global diversity (",
#                                  "estimated to be ", round(divt, 1), " total tips)."),
#                  y="Mean phylogenetic distance to Root", x="Sample size"))
#     plot(slope)
# }; rm(table, clade, data, plateau, div99, frac, divt, slope); dev.off()
# write.table(outData, files$outData, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
# 
