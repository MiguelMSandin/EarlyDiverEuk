#----
#---- loading packages ----

library(ape)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)

#----
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")
rm(list=ls()[!ls() %in% c("files")])
# .rs.restartR()
#----
#---- Diversity fractions --------------------------------------------------------------------------

(files = list(data="plots/fractions/fractions_GoF.tsv",
              outPlotsPref="plots/fractions/fractions"))

if(!dir.exists(dirname(files$outPlot))){dir.create(dirname(files$outPlot))}

data = fread(files$data)

data$alignment = ifelse(grepl("rRA", data$file), "reverse", "forward")
data$root = sub("\\/.*", "", data$file)
data$replicate = data$file %>% sub("\\/clades.*", "", .) %>% sub("root.\\/*", "", .)
data = subset(data, root=="rootA")
data$root = NULL

# data$clade = factor(data$clade, levels=c("Discoba", "Metamonada", 
# 										 "Amoebozoa", "Nucletmycea", "Holozoa", 
# 										 "Haptista", "Cryptista", "Archaeplastida", 
# 										 "Alveolata", "Rhizaria", "Stramenopila"))

datas = data %>% group_by(clade) %>% summarise(meanTrunc=mean(truncated), meanML2l=mean(ml2log2))
datas$mean = (datas$meanTrunc + datas$meanML2l) / 2

datal = data.frame(clade=c(data$clade, data$clade),
                   method=c(rep("truncated", nrow(data)), rep("ML2log", nrow(data))),
                   frac=c(data$truncated, data$ml2log2),
                   tips=c(data$tips, data$tips),
                   tipse=c(data$truncatedTips, data$ml2log2Tips),
                   gofstat=c(data$truncatedGoFstatistic, data$ml2log2GoFstatistic),
                   gofpval=c(data$truncatedGoFpvalue, data$ml2log2GoFpvalue),
                   qqslope=c(data$truncatedQQslope, data$ml2log2QQslope),
                   qqr2=c(data$truncatedQQadjR2, data$ml2log2QQadjR2))

datal$clade = factor(datal$clade, levels=datas$clade[order(datas$mean, decreasing=TRUE)])

(diversities = ggplot(datal, aes(x=clade, y=frac, colour=method))+
    geom_boxplot()+
    geom_jitter(alpha=0.2)+
    # geom_violin(aes(x=clade, y=frac, colour=method), alpha=0.4, scale="width")+
    scale_colour_manual(values=c("#648fff", "#ffb000"))+
    theme_bw()+
    # theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Summary of diversity estimates",
         y="Diversity fraction estimate (percentage of total)", x="Clades"))

(diversitiesGoF = ggplot(datal, aes(x=clade, y=gofstat, colour=method))+
    geom_boxplot()+
    geom_jitter(alpha=0.2)+
    # geom_hline(yintercept=0.05, colour="darkred")+
    # geom_hline(yintercept=0.01, colour="darkred")+
    # geom_violin(aes(x=clade, y=frac, colour=method), alpha=0.4, scale="width")+
    scale_colour_manual(values=c("#648fff", "#ffb000"))+
    theme_bw()+
    # theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Goodness-of-fit Kolmogorov-Smirnov (statistic) distance",
         y="Statistic", x="Clades"))

(diversitiesQQslope = ggplot(datal, aes(x=clade, y=qqslope, colour=method))+
    geom_boxplot()+
    geom_jitter(alpha=0.2)+
    # geom_violin(aes(x=clade, y=frac, colour=method), alpha=0.4, scale="width")+
    scale_colour_manual(values=c("#648fff", "#ffb000"))+
    theme_bw()+
    # theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Q-Q slope",
         y="Q-Q slope", x="Clades"))

(diversitiesQQr2 = ggplot(datal, aes(x=clade, y=qqr2, colour=method))+
    geom_boxplot()+
    geom_jitter(alpha=0.2)+
    # geom_violin(aes(x=clade, y=frac, colour=method), alpha=0.4, scale="width")+
    scale_colour_manual(values=c("#648fff", "#ffb000"))+
    theme_bw()+
    # theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1))+
    labs(title="Q-Q R²",
         y="Q-Q R²", x="Clades"))


pdf(paste0(files$outPlotsPref, ".pdf"), width=11.69, height=8.27*0.6, paper='special'); plot(diversities); dev.off()
pdf(paste0(files$outPlotsPref, "_GoF_statistic.pdf"), width=11.69, height=8.27*0.4, paper='special'); plot(diversitiesGoF); dev.off()
pdf(paste0(files$outPlotsPref, "_qq_slope.pdf"), width=11.69, height=8.27*0.4, paper='special'); plot(diversitiesQQslope); dev.off()
pdf(paste0(files$outPlotsPref, "_qq_r2.pdf"), width=11.69, height=8.27*0.4, paper='special'); plot(diversitiesQQr2); dev.off()


#----
#---- Get rough estimates of eukaryotic diversity --------------------------------------------------

data$estimate = NA
for(n in 1:nrow(data)){
	tmp = subset(data, file==data$file[n] & clade==data$clade[n])
	data$frac[n]
	if(data$frac[n] == min(tmp$frac)){
		data$estimate[n] = "min"
	}else if(data$frac[n] == max(tmp$frac)){
		data$estimate[n] = "max"
	}
}; rm(n, tmp)

data %>% group_by(method) %>% summarise(mean=mean(frac), sd=sd(frac), wmean=weighted.mean(frac, tips))
data %>% group_by(estimate) %>% summarise(mean=mean(frac), sd=sd(frac), wmean=weighted.mean(frac, tips))
data %>% group_by(clade) %>% summarise(mean=mean(frac), sd=sd(frac))

subset(data, estimate=="min") %>% group_by(clade) %>% summarise(mean=mean(frac), sd=sd(frac))
subset(data, estimate=="max") %>% group_by(clade) %>% summarise(mean=mean(frac), sd=sd(frac))



#----
