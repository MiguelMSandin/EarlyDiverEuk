#---- 
#---- loading packages ----

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

#----  
setwd("")
#---- 
#---- Files ----------------------------------------------------------------------------------------

files = list(info ="resources/info_trees.tsv",
			  otus="resources/otu_number_lineage.tsv",
			  plots="plots/phylo")

for(i in 1:length(unlist(strsplit(files$plots, "/")))){
	dir = paste(unlist(strsplit(files$plots, "/"))[1:i], collapse="/")
	if(!dir.exists(dir)){dir.create(dir)}
};rm(i, dir)

#---- 
#---- read info ------------------------------------------------------------------------------------

info = fread(files$info)

{info$CPU_time = as.numeric(info$CPU_time)
info$tree_logLikelihood = as.numeric(info$tree_logLikelihood)
info$tips_in = as.numeric(info$tips_in)
info$tips_out = as.numeric(info$tips_out)
info$pruned = as.numeric(info$pruned)
info$intruders = as.numeric(info$intruders)
info$removed = as.numeric(info$removed)
info$pct = as.numeric(info$pct)}

otus = fread(files$otus)
otus = melt(otus, id.vars=c("alignment", "tree", "replicate"))
colnames(otus) = c("alignment", "tree", "replicate", "group", "OTUs")

tmp = otus %>% group_by(group) %>% summarise(mean=mean(OTUs), sd=sd(OTUs), difference=max(OTUs)-min(OTUs))
otus = merge(otus, tmp, by="group")
otus$group = factor(otus$group, levels=unique(otus$group[order(otus$mean)]))

summ = otus %>% group_by(group) %>% summarise(difference=max(OTUs)-min(OTUs), relative=(max(OTUs)-min(OTUs))/mean(OTUs)*100)
(summ = summ[order(summ$difference, decreasing=TRUE),])

#---- 
#---- Plot Log-likelihood of all different phylogenetic analyses -----------------------------------

infof = subset(info, !is.na(info$pruned))
ggplot(infof, aes(x=step, y=tree_logLikelihood, colour=file))+
	geom_jitter(alpha=0.6, height=0)+
	geom_boxplot(alpha=0.4)+
	theme_bw()

(likeLsteps = ggplot(infof, aes(x=step, y=tree_logLikelihood, colour=file))+
        geom_jitter(alpha=0.6, height=0)+
        geom_boxplot(alpha=0.4)+
        facet_wrap(~step, scales="free")+
        theme_bw())

pdf(paste0(files$plots, "/likelihood_vs_aligment.pdf"), width=11.69, height=8.27, paper='special'); plot(likeLsteps); dev.off()

#---- 
#---- OTU number per lineage -----------------------------------------------------------------------

ggplot(otus, aes(x=group, y=log(OTUs)))+
	geom_point()+
	theme_bw()+
	theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

(boxes = ggplot(otus, aes(x=group, y=log(OTUs)))+
		geom_jitter(alpha=0.6, position=position_jitter(0.02))+
		geom_boxplot(alpha=0.4)+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))

pdf(paste0(files$plots, "/OTUnumber_per_group_variability.pdf"), width=11.69, height=8.27, paper='special'); plot(boxes); dev.off()

(boxWrap = ggplot(otus, aes(x=group, y=OTUs))+
		geom_jitter(alpha=0.6, position=position_jitter(0.02))+
		geom_boxplot(alpha=0.4)+
		facet_wrap(~group, scales="free")+
		theme_bw())

pdf(paste0(files$plots, "/OTUnumber_per_group_variability_wrapped.pdf"), width=11.69, height=8.27, paper='special'); plot(boxWrap); dev.off()

# And now relative to the total number of reads in the tree
tmp = otus %>% group_by(group) %>% mutate(relToMax=OTUs/max(OTUs))
tmp$group = factor(tmp$group, levels=unique(tmp$group[order(tmp$relToMax)]))
tmps = tmp %>% select(c("group", "mean", "sd"))
tmps = tmps[!duplicated(tmps),]

(relToMaxPlot = ggplot(tmp, aes(x=group, y=relToMax))+
		geom_jitter(height=0)+
		geom_boxplot(alpha=0.4)+
		geom_text(data=tmps, aes(label=paste0("Average # tips\n", round(mean), " (sd:", round(sd, 1), ")"), y = -0.05), 
				  size=3, position=position_dodge(width=0.8))+
		theme_bw()+
		theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))

pdf(paste0(files$plots, "/OTUnumber_per_group_variability_relativeToTotal.pdf"), width=11.69, height=8.27, paper='special'); plot(relToMaxPlot); dev.off()

#---- 
