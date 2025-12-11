#----
#---- loading packages and data.frames -------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(treeio)
library(tidytree)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----
#---- set working directory and file names ---------------------------------------------------------
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")

rm(list=ls()[!ls() %in% c("geo", "files")])
# .rs.restartR()

files = list(extTrees="LTT-data.tsv",
             dirClades="clades/",
             extClades="Archaeplastida_annot_LTT.tsv",
             factor=-1,
             out="plots/LTT/MC01_LTT_Archaeplastida") 

#----
#---- Read tables ----------------------------------------------------------------------------------

# First search all the tables in the directories recursively
files$tables = grep(files$extTrees, dir(recursive=TRUE), value=TRUE)
files

# Concatenate all tables into one
data = data.frame()
i = 0
for(file in files$tables){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tables), sep="")
	f = fread(paste0(file))
	f = subset(f, tree=="Main_tree" | tree=="Archaeplastida" | tree=="Discoba" | tree=="Amoebozoa" | tree=="Rhizaria" | tree=="Alveolata" | tree=="Nucletmycea" | tree=="Holozoa")
	f$root = gsub("/.*", "", file)
	f$run = file %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)
	f$hpd_max = NULL
	f$hpd_min = NULL
	data = rbind(data, f)
}; rm(file, f, i); cat("\n")

# Now read the LTT plots of the Archaeplastida only
files$clades = grep(files$extClades, dir(recursive=TRUE), value=TRUE)
files$clades
i = 0
for(file in files$clades){
  cat("\r  Reading files ", (i = i + 1), "/", length(files$clades), sep="")
  f = fread(paste0(file))
  f$time = -f$time
  colnames(f) = c("tree", "time", "N", "logN")
  f = subset(f, tree=="Chloroplastida" | tree=="Rhodophyta")
  f$root = gsub("/.*", "", file)
  f$run = file %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)
  data = rbind(data, f)
}; rm(file, f, i); cat("\n")
data$tree = sub("main", "Main_tree", data$tree)

# # Or the trees instead
# files$clades = grep(files$extClades, dir(recursive=TRUE), value=TRUE)
# files$clades
# 
# i = 0
# for(file in files$clades){
# 	cat("\r  Reading tree\t", (i = i + 1), "/", length(files$clades), sep="")
# 	treei = read.beast(file)
# 	treed = as_tibble(treei)
# 	names(treed) = gsub("!", "", names(treed))
# 	annotations = sort(as.vector(treed$name[!is.na(treed$name)]))
# 	for(annot in annotations){
# 		node = as.numeric(subset(treed, name==annot)$node)
# 		tmp = as.phylo(tree_subset(treei, node=node, levels_back=0))
# 		if(Ntip(tmp) > 300){
# 			tmp = as.data.frame(ltt.plot.coords(tmp))
# 			tmp = data.frame(tree=annot,
# 							 time=as.numeric(tmp$time),
# 							 N=tmp$N,
# 							 logN=log(tmp$N, exp(1)),
# 							 root=gsub("/.*", "", file),
# 							 run=file %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .))
# 			data = rbind(data, tmp)
# 		}
# 	}
# }; rm(file, i, treei, treed, tmp, annot, node); cat("\n")

write.table(data, paste0(files$out, ".tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#----
#---- Summarise the table --------------------------------------------------------------------------

tmp = data %>% group_by(root, tree) %>% reframe(taxa=length(unique(run)))
all(tmp$taxa==32)

summ = subset(data, tree!="Chlorophyta" & tree!="Streptophyta") %>% 
  group_by(tree, root, N, logN) %>% 
  summarise(median=quantile(time, 0.50),
            q05=quantile(time, 0.05),
            q95=quantile(time, 0.95),
            replicates=length(unique(paste(root, run))))

write.table(summ, paste0(files$out, "_summary.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#----
#---- beautifying the dataset for plotting ---------------------------------------------------------

sort(unique(summ$tree))

summ$tree = factor(summ$tree, 
				   levels=c("Main_tree", "Archaeplastida", "Rhodophyta", "Chloroplastida", # "Chlorophyta", "Streptophyta",
				   		 "Discoba", "Rhizaria", "Alveolata", "Amoebozoa", "Nucletmycea", "Holozoa"))

{summ$colour = as.character(summ$tree)
	summ$colour[which(summ$colour=="Main_tree")]="black"
	summ$colour[which(summ$colour=="Archaeplastida")]="darkseagreen3"
	summ$colour[which(summ$colour=="Rhodophyta")]="red4"
	summ$colour[which(summ$colour=="Chloroplastida")]="olivedrab"
	summ$colour[which(summ$colour=="Chlorophyta")]="olivedrab3"
	summ$colour[which(summ$colour=="Streptophyta")]="olivedrab1"
	summ$colour[which(summ$colour=="Discoba")]="orange2"
	summ$colour[which(summ$colour=="Rhizaria")]="darkorchid2"
	summ$colour[which(summ$colour=="Alveolata")]="darkorchid4"
	summ$colour[which(summ$colour=="Amoebozoa")]="royalblue1"
	summ$colour[which(summ$colour=="Nucletmycea")]="steelblue2"
	summ$colour[which(summ$colour=="Holozoa")]="steelblue4"}

summ$colour = factor(summ$colour, 
				   levels=c("black", "darkseagreen3", "red4", "olivedrab", "olivedrab3", "olivedrab1",
				   		 "orange2", "darkorchid2", "darkorchid4", "royalblue1", "steelblue2", "steelblue4"))

summ$root = factor(summ$root, levels=c("rootD", "rootA"))

#----
#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table

geos = subset(geo, time > min(data$time, na.rm=TRUE))

(lttplot = ggplot(summ)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_ribbon(aes(xmin=q95, xmax=q05, y=N, fill=tree), alpha=0.2, colour = NA) +
		geom_line(aes(x=median, y=N, color=tree)) +
		geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
		scale_x_continuous(breaks=seq((round(min(summ$q95, na.rm=TRUE), -2)), 0, 100), 
						   minor_breaks=seq((round(min(summ$q95, na.rm=TRUE), -2)), 0, 50)) +
		scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid( ~ root, scales="free")+
		scale_color_manual(values=as.character(sort(unique(summ$colour))))+
		scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
		theme_classic()+
		labs(y="Ln (N)", x="Time (Ma)"))

pdf(paste0(files$out, ".pdf"), width=11.69, height=8.27, paper='special'); plot(lttplot); dev.off()

#----
