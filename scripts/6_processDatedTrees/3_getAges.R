#---- 
#---- loading packages ----

library(ape)
library(treeio)

library(data.table)
library(tidyr)
library(dplyr)

library(ggplot2)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				 era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----  
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/")
#---- 
#---- set file names and directories ---------------------------------------------------------------

rm(list=ls()[!ls() %in% c("files")])

files = list(extension="median\\.tre",
             # extension="median_LTT-data\\.tsv",
             out="plots/ages/",
             outPrefix="MC01-MC02")

if(!dir.exists(dirname(files$out))){dir.create(dirname(files$out))}
#---- 
#---- 
#----  Get ages of all annotated selected trees ----------------------------------------------------

files$trees = grep(paste0(files$extension, "$"), dir(recursive=TRUE), value=TRUE)
(files$trees = grep("MC0(1|2)", files$trees, value=TRUE))

ages = data.frame()
for(tree in files$trees){
    cat("Reading tree '", tree, "' (", grep(tree, files$trees), "/", length(files$trees), ")\n", sep="")
    treei = read.beast(tree)
    # Transform files
    treed = as_tibble(treei)
    names(treed) = gsub("!", "", names(treed))
    # Select all annotated nodes
    annotations = sort(as.vector(treed$name[!is.na(treed$name)]))
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
        }else if (annot=="Alveolata" | annot=="Amoebozoa" | annot=="Archaeplastida" |
        		  annot=="Cryptista" | annot=="Discoba" | annot=="Haptista" | annot=="Holozoa" |
        		  annot=="Metamonada" | annot=="Nucletmycea" | annot=="Rhizaria" | annot=="Stramenopila"){
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
# ages$tree = sub("^stepDating\\/", "", ages$tree)
ages$calibration = ages$tree %>% sub("\\/root.*", "", .) %>% sub("-2", "", .)
ages$root = ages$tree %>% sub("MC0(1|2)-2\\/","",.) %>% sub("\\/(r|f)RA.*","",.)
ages$alignment = ifelse(grepl("rRA", ages$tree), "reverse", "forward")
ages$replicate = ages$tree %>% sub(".*step3._RA","",.) %>% sub("_rootA.*|_rootD.*", "", .) %>% sub("_rep", "-", .)

ages$height = as.numeric(ages$height)

table(ages$clade)
table(ages$replicate)
#

write.table(ages, paste0(files$out, "/", files$outPrefix,"_clade_ages.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#----
#----  Or get ages of all exported LTT tables from the trees ---------------------------------------

# files$ltts = grep(paste0(files$extension, "$"), dir(recursive=TRUE), value=TRUE)
# 
# ages = data.frame()
# for(ltt in files$ltts){
#   cat("Reading LTT data '", ltt, "' (", grep(ltt, files$ltts), "/", length(files$ltts), ")\n", sep="")
#   lttd = fread(ltt)
#   lttd = subset(lttd, N==1)
#   for(annot in unique(lttd$tree)){
#     cat("\r    Working on ", annot, " (", grep(annot, unique(lttd$tree)), "/", length(unique(lttd$tree)), ")                   ", sep="")
#     if(annot=="Main_tree"){
#       tmp = treed[is.na(treed$label),][1,]
#       ages = rbind(ages, data.frame(tree=tree,
#                                     clade="main",
#                                     height=tmp$height_median,
#                                     HPDmin=tmp$height_0.95_HPD[[1]][1],
#                                     HPDmax=tmp$height_0.95_HPD[[1]][2]))
#     }else if (annot=="Alveolata" | annot=="Amoebozoa" | annot=="Archaeplastida" | 
#               annot=="Cryptista" | annot=="Discoba" | annot=="Haptista" | annot=="Holozoa" | 
#               annot=="Metamonada" | annot=="Nucletmycea" | annot=="Rhizaria" | annot=="Stramenopila"){
#       tmp = subset(treed, name==annot)
#       ages = rbind(ages, data.frame(tree=tree,
#                                     clade=annot,
#                                     height=tmp$height_median,
#                                     HPDmin=tmp$height_0.95_HPD[[1]][1],
#                                     HPDmax=tmp$height_0.95_HPD[[1]][2]))
#     }
#   }
#   rm(annot, tmp)
#   cat("\r    Done                                   \n")
# };rm(tree)
# 
# unique(ages$tree)
# # ages$tree = sub("^stepDating\\/", "", ages$tree)
# ages$alignment = ifelse(grepl("rRA", ages$tree), "reverse", "forward")
# ages$root = ages$tree %>% sub("MC\\d+\\/","",.) %>% sub("MC\\d+r\\/","",.) %>% sub("\\/.*", "", .)
# ages$replicate = ages$tree %>% sub(".*step3._RA","",.) %>% sub("_rootA.*|_rootD.*", "", .) %>% sub("_rep", "-", .)
# 
# ages$height = as.numeric(ages$height)
# 
# table(ages$clade)
# table(ages$replicate)
# #
# 
# write.table(ages, paste0(files$out,"/clade_ages.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#----
#---- Read table -----------------------------------------------------------------------------------

ages = fread(paste0(files$out, "/", files$outPrefix,"_clade_ages.tsv"))

#----
#---- Subsetting -----------------------------------------------------------------------------------

ages = fread(paste0(files$out,"_clade_ages.tsv"))

ages = subset(ages, clade!="Ancyromonadida" & clade!="Apusomonadida"  & clade!="Apusomonaida" & clade!="Breviatea" & 
                   clade!="CRuMs" & clade!="CRUMs" & clade!="Hemimastigophora" & clade!="Malawimonadida" & clade!="Telonemia" & !is.na(clade))
ages$clade = factor(ages$clade, levels=c("main", 
                                         "Discoba", "Metamonada", 
                                         "Amoebozoa", "Nucletmycea", "Holozoa", 
                                         "Haptista", "Cryptista", "Archaeplastida", 
                                         "Alveolata", "Rhizaria", "Stramenopila"))

summ = ages %>% group_by(clade, calibration, root) %>% summarise(median=median(height), #HPDmin=mean(HPDmin), HPDmax=mean(HPDmax),
													 max=max(height), min=min(height))

(summ[order(summ$median, decreasing=TRUE),])

#----
#---- Plotting -------------------------------------------------------------------------------------

tmp = unique(agess$clade)

pdf(paste0(files$out, "/", files$outPrefix,"_clade_ages.pdf"), width=11.69, height=8.27, paper='special')
for(node in unique(ages$clade)){
    cat("\r  Plotting ", node, " (", grep(node, unique(ages$clade)), "/", length(unique(ages$clade)), ")                    ", sep="")
    tmp = subset(ages, clade==node)
    
    agesPlot = ggplot(tmp, aes(x=replicate, y=height, colour=root))+
      geom_segment(aes(x=replicate, xend=replicate, y=HPDmin, yend=HPDmax, colour=root), alpha=0.4, linewidth=2, lineend="round")+
      geom_point() +
      # facet_wrap(~calibration+root+alignment, nrow=1, scales="free_x")+
      facet_grid(~calibration+root+alignment, scales="free_x")+
      # facet_wrap(~clade, scales="free")+
      scale_y_continuous(breaks=seq(0, round(max(ages$HPDmax), -2), 100),
                         limits=c(0, round(max(ages$HPDmax), -2))) +
      theme_classic()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      labs(title=paste0(node),
           y="Time (Ma)", x="Replicates")
    plot(agesPlot)
}; rm(node, tmp, agesPlot); cat("\rDone                                        \n"); dev.off()

#
#---- Plot now a summary ---------------------------------------------------------------------------

unique(ages$clade)
tmp = ages %>% group_by(clade) %>% summarise(median=median(height))
tmp = tmp$clade[order(tmp$median, decreasing=TRUE)]
ages$clade = factor(ages$clade, levels=tmp)

geos = subset(geo, time < max(ages$HPDmax, na.rm=TRUE))

ages$grouping = paste0(ages$calibration, "-", ages$root)
ages$grouping = factor(ages$grouping, levels=c("MC01-rootD", "MC01-rootA", "MC02-rootD", "MC02-rootA"))

(agesSumPlot = ggplot(ages)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_point(aes(x=-HPDmin, y=clade, group=grouping), position = position_dodge(0.4), alpha=0.4, colour="grey80") +
		geom_point(aes(x=-HPDmax, y=clade, group=grouping), position = position_dodge(0.4), alpha=0.4, colour="grey80") +
		geom_point(aes(x=-height, y=clade, colour=grouping), position = position_dodge(0.4), alpha=0.4) +
		geom_violin(aes(x=-height, y=clade, colour=grouping), position = position_dodge(0.4), alpha=0.4) +
		# geom_boxplot(aes(x=-height, y=clade, colour=root), position = position_dodge(0.4), alpha=0.4) +
		# facet_wrap(~root, nrow=1)+
		scale_x_continuous(breaks=seq(round(min(-ages$HPDmax), -2), 0, 100),
						   limits=c(round(min(-ages$HPDmax), -2), 0)) +
		theme_classic()+
    theme(axis.title.y=element_blank())+
    labs(x="Time (Ma)"))

pdf(paste0(files$out, "/", files$outPrefix, "_clade_ages_summary.pdf"), width=11.69, height=8.27, paper='special'); plot(agesSumPlot); dev.off()

#----
