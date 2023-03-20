#----
#---- loading packages ----

library(ape)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)

#----
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/")
#----
#---- Diversity fractions --------------------------------------------------------------------------

(files = list(files=grep("fractions\\/.*fractions\\.tsv", dir(recursive=TRUE), value=TRUE),
			  outPlot="plots/fractions/fractions.pdf"))

if(!dir.exists("plots/fractions")){dir.create("plots/fractions")}

data = data.frame()
for(file in files$files){
	cat("\r  Reading file '", file, "' (", grep(file, files$files), "/", length(files$files), ")", sep="")
	tmp = fread(file)
	tmp$file = file
	data = rbind(data, tmp)
}; rm(file, tmp); cat("\rDone", rep(" ", 100), "\n")

data$alignment = ifelse(grepl("^MC01r", data$file), "reverse", "forward")
data$root = data$file %>% sub("MC\\d+\\/","",.) %>% sub("MC\\d+r\\/","",.) %>% sub("\\/.*", "", .)
data$replicate = data$file %>% sub(".*fractions\\/RA","",.) %>% sub("_supergroups.*", "", .)

data$clade = factor(data$clade, levels=c("Discoba", "Metamonada", 
										 "Amoebozoa", "Nucletmycea", "Holozoa", 
										 "Haptista", "Cryptista", "Archaeplastida", 
										 "Alveolata", "Rhizaria", "Stramenopila"))

(diversities = ggplot(data, aes(x=clade, y=frac, colour=method, size=tips))+
		geom_point(alpha=0.4)+
		theme_bw()+
		theme_minimal()+
		theme(axis.text.x=element_text(angle=30, hjust=1))+
		labs(title="Summary of diversity estimates",
			 y="Diversity fraction estimate (percentage of total)", x="Clades"))

pdf(paste0(files$outPlot), width=11.69, height=8.27, paper='special'); plot(diversities); dev.off()

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
	}else{
		data$estimate[n] = "med"
	}
}; rm(n, tmp)

data %>% group_by(estimate) %>% summarise(mean=mean(frac), sd=sd(frac))
tmp = data %>% group_by(clade, estimate) %>% summarise(mean=mean(frac), sd=sd(frac))

subset(data, estimate=="min") %>% group_by(clade) %>% summarise(mean=mean(frac), sd=sd(frac))
subset(data, estimate=="max") %>% group_by(clade) %>% summarise(mean=mean(frac), sd=sd(frac))

#----
