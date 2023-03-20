#---- 
#---- loading packages ----

library(ape)
library(phangorn)
library(treeio)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)

#----  
setwd("")
rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()
#---- 
#---- Set names of trees, files and output ---------------------------------------------------------

(files = list(
	treeDated=grep("_median\\.tre$", dir(), value=TRUE), # The dated tree
	treePhylo=grep("_iqtreef.*\\.tre$", dir(), value=TRUE), # The phylogenetic tree
	abundances="../../../../resources/tipNames_OTUreads_norm.tsv", # The abundances table
	dirClades="clades", # The output directory where all subclades of the dated tree are going to be exported
	dirTrees="clades/phylo", # The output directory where all subclades of the phylogenetic tree are going to be exported	
	outDir="fractions",  # The output directory where all diversity fractions are going to be exported
	out=paste0(sub(".*\\/", "", getwd()), "_supergroups"), # A prefix for the output
	heightFactor=100, # A factor to multiply the branch-lengths of the dated tree if needed
	minTips=300, # A minimum number of tips to consider or not a supergroup
	maxSamplingFactor=10)) # The maximum number of rarefied samples (maxSamplingFactor x number of tips)

if(!dir.exists(files$dirClades)){dir.create(files$dirClades)}
if(!dir.exists(files$dirTrees)){dir.create(files$dirTrees)}
if(!dir.exists(files$outDir)){dir.create(files$outDir)}

#---- 
#---- Transform scale of dated tree if needed and export independent clades ------------------------

# Read annotated nexus tree file
treei <- read.beast(files$treeDated)

# Transform tree into a data.frame
treed <- as_tibble(treei)
names(treed) <- gsub("!", "", names(treed))

hist(treed$branch.length)
# hist(log(treed$branch.length))
summary(treed$branch.length)

# Transform scale if needed
# tree = read.tree(file)
# tree$edge.length = tree$edge.length * files$heightFactor
# write.tree(tree, sub("\\.[^\\.]+$", "_timeScaled.tre", file))
# treed <- as_tibble(treei)
# names(treed) <- gsub("!", "", names(treed))

# Select all annotated nodes
(annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))

# Check for misspellings of the 11 most abundant clades in the annotations
tmp = c("Discoba", "Metamonada", "Amoebozoa", "Nucletmycea", "Holozoa", "Haptista", "Cryptista", "Archaeplastida", "Rhizaria", "Stramenopila", "Alveolata")
if(!all(tmp %in% annotations)){tmp[!tmp %in% annotations]}else{cat("All major clades are correctly spelled in the tree")}

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

#---- Export independent clades from phylogenetic tree --------------------------------------------

rm(list=ls()[!ls() %in% c("files")])

# Read annotated nexus tree file
treei <- read.beast(files$treePhylo)

# Transform tree into a data.frame
treed <- as_tibble(treei)
names(treed) <- gsub("!", "", names(treed))

hist(treed$branch.length)
# hist(log(treed$branch.length))
summary(treed$branch.length)

# Select all annotated nodes
(annotations <- sort(as.vector(treed$name[!is.na(treed$name)])))

# Check for misspellings of the 11 most abundant clades in the annotations
tmp = c("Discoba", "Metamonada", "Amoebozoa", "Nucletmycea", "Holozoa", "Haptista", "Cryptista", "Archaeplastida", "Rhizaria", "Stramenopila", "Alveolata")
if(!all(tmp %in% annotations)){tmp[!tmp %in% annotations]}else{cat("All major clades are correctly spelled in the tree")}

# loop through all annotated clades, and export those with at least 'files$minTips' tips
for(annot in annotations){
    cat("  Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")\n", sep="")
    node = as.numeric(subset(treed, name==annot)$node)
    tmp = tree_subset(treei, node=node, levels_back=0)
    tmp = as.phylo(tmp)
    if(Ntip(tmp) >= files$minTips){
        cat("    Exporting tree with", Ntip(tmp),"tips\n")
        write.tree(tmp, paste0(files$dirTrees, "/clade_", annot, ".tre"))
    }else{cat("    ! Clade not exported (", Ntip(tmp), " tips)\n", sep="")}
};rm(annot, node, tmp)

#

#---- Check all exported clades from phylogenetic and dated tree have the same Ntip ----------------

rm(list=ls()[!ls() %in% c("files")])

tmp1 = grep("^clade_.*\\.tre$", dir(files$dirClades), value=TRUE)
tmp2 = grep("^clade_.*\\.tre$", dir(files$dirTrees), value=TRUE)

# First check we have all clades in both the phylogenetic and dated tree
all(tmp1 %in% tmp2)

# Now check the number of tips
for(tree in tmp1){
	clade = tree %>% sub("clade_", "", .) %>% sub("\\..*", "", .)
	phylo = read.tree(paste0(files$dirTrees, "/", tree))
	dated = read.tree(paste0(files$dirClades, "/", tree))
	if(	Ntip(phylo) == Ntip(dated) ){
		cat("  Clade", clade, "has", Ntip(phylo), "in both trees: \tOK\n")
	}else{
		cat("!!Clade", clade, "has", Ntip(dated), "in dated tree and", Ntip(phylo), "in phylogenetic tree...\n")
	}
}; rm(tree, clade, phylo, dated, tmp1, tmp2)

#

#----
#---- Estimate phylogenetic species richness based on the rarefied trees ---------------------------

rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()

# Read abundances file
abun = fread(files$abundances)

# Rarefy the trees _________________________________________________________________________________
pdf(paste0(files$outDir, "/", files$out, "_rarefiedPhylo.pdf"), width=11.69, height=8.27, paper='special')
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
        	tips = sample(dists$tips, i, replace=TRUE, prob=dists$abunRel)
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
    
    outData = rbind(outData, data.frame(clade= clade, 
                                        tips=Ntip(phylo), 
                                        fraction=frac,
                                        total=divt))
    
    # plot
    slope = ggplot(rares, aes(x=sample, y=mean))+
        geom_line(aes(y=min), colour="lightgrey")+
        geom_line(aes(y=max), colour="lightgrey")+
        geom_hline(yintercept=total, color="springgreen3")+
        geom_segment(aes(x=Ntip(phylo), y=0, 
                         xend=Ntip(phylo), yend=div1), colour="springgreen3")+
        geom_segment(aes(x=0, y=div1, 
                         xend=Ntip(phylo), yend=div1), colour="springgreen3")+
        geom_line()+
        theme_minimal()+
        labs(title=paste0(clade, ": ", Ntip(phylo), " tips"),
             subtitle=paste0("This tree represents an estimated ~", round(frac), "% of the global diversity (",
                             "estimated to be ~", round(divt), " total tips) based on the total phylogenetic diversity.\n"),
             y="Mean phylogenetic distance to Root", x="Sample size")
    plot(slope)
}; rm(tree, phylo, rare, rares, dists, total, abuns, div1, frac, divt, clade, k, samples, i, j, tmp, slope); dev.off()

write.table(outData, paste0(files$outDir, "/", files$out, "_rarefiedPhylo.tsv"), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(outSlopes, paste0(files$outDir, "/", files$out, "_rarefiedPhylo_raw.tsv"), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#

#---- Estimate species richness based on reads abundances ------------------------------------------

rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()

abundances = fread(files$abundances)

trees = grep("^clade_.*\\.tre$", dir(paste0(files$dirTrees)), value=TRUE); i = 0
dataOut = data.frame()
pdf(paste0(files$outDir, "/", files$out, "_preston.pdf"), width=11.69, height=8.27, paper='special')
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

write.table(dataOut, paste0(files$outDir, "/", files$out, "_preston.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
#---- Check and compare the different estimates ----------------------------------------------------

rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()

tmp1 = fread(paste0(files$outDir, "/", files$out, "_preston.tsv"))
tmp2 = fread(paste0(files$outDir, "/", files$out, "_rarefiedPhylo.tsv"))

if(all(tmp1$clade == tmp2$clade) & all(tmp1$tips == tmp2$tips)){
    data = data.frame(clade=rep(tmp1$clade, 3),
                      tips=rep(tmp1$tips, 3),
                      method=c(rep("quasiPoisson", nrow(tmp1)), rep("maxLikelihood2log2", nrow(tmp1)), rep("rarefied", nrow(tmp2))),
                      frac=c(tmp1$quasiPoisson, tmp1$maxLikelihood2log2, tmp2$fraction),
                      total=c(tmp1$quasiPoissonTips, tmp1$maxLikelihood2log2Tips, tmp2$total))
    rm(tmp1, tmp2)
}else{cat("\nPlease check that the clades and tips in both files are in matching order\n\n")}

write.table(data, paste0(files$outDir, "/", files$out, "_fractions.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

unique(data$clade)
data$clade = factor(data$clade, levels=c("Discoba", "Metamonada","Amoebozoa", "Nucletmycea", "Holozoa","Haptista",  "Cryptista", "Archaeplastida","Rhizaria", "Stramenopila","Alveolata"))

(summ = ggplot(data, aes(x=clade, y=frac))+
    geom_point(aes(colour=method, size=tips))+
    theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Summary of diversity estimates",
         y="Diversity fraction estimate (percentage of total)", x="Clades"))

pdf(paste0(files$outDir, "/", files$out, "_fractions.pdf"), width=11.69, height=8.27, paper='special'); plot(summ); dev.off()

#----
