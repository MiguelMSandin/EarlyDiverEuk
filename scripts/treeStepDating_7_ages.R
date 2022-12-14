#---- 
#---- loading packages ----

library(ape)
library(treeio)

library(BAMMtools)
library(HDInterval)

library(data.table)
library(tidyr)
library(dplyr)

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
setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/")
#---- 
#---- set file names and directories ---------------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirs=c("MC01/rootA/", "MC01/rootD/", "MC01r/rootA/", "MC01r/rootD/"),
             out="plots/ages")

#---- 
#---- Read trees  ----------------------------------------------------------------------------------

files$trees = c()
for(dir in files$dirs){
    files$trees = c(files$trees, paste0(dir, grep("median\\.tre$", dir(dir, recursive=TRUE), value=TRUE)))
}; rm(dir)


#---- 
#----  Get ages of all annotated selected trees ----------------------------------------------------

ages = data.frame()
for(tree in files$trees){
    cat("Reading tree '", tree, "' (", grep(tree, files$trees), "/", length(files$trees), ")\n", sep="")
    treei <- read.beast(tree)
    
    # Transform files
    treed <- as_tibble(treei)
    names(treed) <- gsub("!", "", names(treed))
    
    # Select all annotated nodes
    annotations <- sort(as.vector(treed$name[!is.na(treed$name)]))
    
    # Loop through the main tree and the different annotated nodes to get the ages
    cat("  Extracting ages\n")
    for(annot in c("Main_tree", annotations)){
        cat("\r    Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")                   ", sep="")
        if(annot=="Main_tree"){
            tmp = treed[is.na(treed$label),][1,]
            ages = rbind(ages, data.frame(tree=tree,
                                          clade="main",
                                          height=tmp$height_median,
                                          HPDmin=tmp$height_0.95_HPD[[1]][1],
                                          HPDmax=tmp$height_0.95_HPD[[1]][2]))
        }else{
            tmp = subset(treed, name==annot)
            ages = rbind(ages, data.frame(tree=tree,
                                          clade=annot,
                                          height=tmp$height_median,
                                          HPDmin=tmp$height_0.95_HPD[[1]][1],
                                          HPDmax=tmp$height_0.95_HPD[[1]][2]))
        }
    }
    rm(annot, tmp)
    cat("\r    Done                                   \n")
};rm(tree)

unique(ages$tree)
ages$alignment = ifelse(grepl("^MC01r", ages$tree), "reverse", "forward")
ages$root = ages$tree %>% sub("MC\\d+\\/","",.) %>% sub("MC\\d+r\\/","",.) %>% sub("\\/.*", "", .)
ages$replicate = ages$tree %>% sub(".*step3._RA","",.) %>% sub("_rootA.*|_rootD.*", "", .) %>% sub("_rep", "-", .)

ages$height = as.numeric(ages$height)
#

write.table(ages, paste0(files$out,"/clade_ages.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

#----
#---- Subsetting -----------------------------------------------------------------------------------

agess = subset(ages, clade!="Ancyromonadida" & clade!="Apusomonadida"  & clade!="Apusomonaida" & clade!="Breviatea" & 
                   clade!="CRuMs" & clade!="CRUMs" & clade!="Hemimastigophora" & clade!="Malawimonadida" & clade!="Telonemia" & !is.na(clade))
agess$clade = factor(agess$clade, levels=c("main", 
                                         "Discoba", "Metamonada", 
                                         "Amoebozoa", "Nucletmycea", "Holozoa", 
                                         "Haptista", "Cryptista", "Archaeplastida", 
                                         "Alveolata", "Rhizaria", "Stramenopila"))

summ = agess %>% group_by(clade, alignment, root) %>% summarise(mean=mean(height), HPDmin=mean(HPDmin), HPDmax=mean(HPDmax))

#----
#---- Plotting -------------------------------------------------------------------------------------

pdf(paste0(files$out,"/clade_ages.pdf"), width=11.69, height=8.27, paper='special')
for(node in unique(agess$clade)){
    cat("\r  Plotting '", node, "' (", grep(node, unique(agess$clade)), "/", length(unique(agess$clade)), ")                    ", sep="")
    tmp = subset(agess, clade==node)
    tmps = subset(summ, clade==node)
    agesPlot <- ggplot(tmp, aes(x=replicate, y=height, colour=root))+
        geom_segment(aes(x=replicate, xend=replicate, y=HPDmin, yend=HPDmax, colour=root), alpha=0.4, size=2, lineend="round")+
        geom_point() +
        facet_wrap(~alignment+root, nrow=1, scales="free_x")+
        # facet_wrap(~clade, scales="free")+
        scale_y_continuous(breaks=seq(0, round(max(agess$HPDmax), -2), 100),
                           limits=c(0, round(max(agess$HPDmax), -2))) +
        theme_classic()+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        labs(title=paste0(node, "; median: ", round(median(tmp$height)), " (min: ", round(min(tmp$height)), "- max: ", round(max(tmp$height)), "; SDmedian: ", round(sd(tmp$height)),") Ma"),
             subtitle=paste0("HPDmin: ", round(min(tmp$HPDmin)), "; HPDminMedian: ", round(median(tmp$HPDmin)), " - HPDmaxMedian: ", round(median(tmp$HPDmax)), "; HPDmax: ", round(max(tmp$HPDmax)), " Ma"),
             y="Time (Ma)", x="Replicates")
    agesPlot
    plot(agesPlot)
}; rm(node, tmp, tmps, agesPlot); cat("\rDone                                        "); dev.off()

#
#----
#----