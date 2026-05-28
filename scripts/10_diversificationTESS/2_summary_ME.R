#----
#---- Loading packages -----------------------------------------------------------------------------

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

geo = data.frame(time= c(2500, 2300, 2050, 1800, 1600, 1400, 1200, 1000, 720, 635, 541, 485.4, 443.8, 419, 358.9, 298.9, 251.9, 201.4, 145, 66, 23.03, 2.58),
                 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                 era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

era = data.frame(time= c(2500, 1600, 1000, 541, 251.9, 66),
                 era=c("Paleoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Paleozoic", "Mesozoic", "Cenozoic"))

#----
#---- Setting working directory and file names -----------------------------------------------------

setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")
rm(list=ls()[!ls() %in% c("geo", "era", "files")])
# .rs.restartR()

files = list(files=grep("comet.*_d.*.tsv", dir(recursive=TRUE), value=TRUE),
             outPrefix="plots/massExtinctions/mass_extinctions_CoMET")

(files$root = unique(sub("\\/.*", "", files$files)))
(files$replicates = unique(files$files %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)))
(files$clades = unique(files$files %>% sub(".*\\/", "", .) %>% sub("_d.*", "", .)))

for(i in files$root){cat("  ", i, "\t", sum(grepl(i, files$files)), "\n")}
for(i in files$replicates){cat("  ", i, "\t", sum(grepl(i, files$files)), "\n")}
for(i in files$clades){cat("  ", i, "\t", sum(grepl(i, files$files)), "\n")}
rm(i)

#----
#---- Reading files --------------------------------------------------------------------------------

# Read files
data = data.frame()
for(file in files$files){
  cat("\r  Reading file ", grep(file, files$files), "/", length(files$files), sep="")
  tmp = fread(file)
  tmp$episode = 1:nrow(tmp)
  root = sub("\\/.*", "", file)
  replicate = file %>% sub("root.\\/", "", .) %>% sub("\\/comet.*", "", .)
  clade = file %>% sub(".*comet\\/", "", .) %>% sub("\\/.*", "", .)
  diversity = as.numeric(file %>% sub(".*_d", "", .) %>% sub("\\..*", "", .))
  tmp = cbind(root, replicate, clade, diversity, tmp)
  data = rbind(data, tmp)
}; rm(file, tmp, root, replicate, clade, diversity); cat("\n")

# Categorizing diversity fractions to "minimum" and "maximum" scenarios
data$div = NA
for(cladei in unique(data$clade)){
  div = subset(data, clade==cladei)$diversity
  m = mean(div)
  data$div[which(data$clade==cladei & data$diversity > m)] = "max"
  data$div[which(data$clade==cladei & data$diversity < m)] = "min"
};rm(cladei, div, m)
if(any(is.na(data$div))){cat("  Warning! There are some diversity fractions not categorized. Please check!")}

#----
#---- Exporting data -------------------------------------------------------------------------------

write.table(data, paste0(files$outPrefix, ".tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#----
#---- Summarizing tables individually per clades ---------------------------------------------------

# data$bfp = ifelse(data$bf >= 0, data$bf, NA) 
# datas = data %>% group_by(root, clade, div, episode) %>% summarise(bfmax=max(bf), 
#                                                                    bfmean=mean(bf),
#                                                                    bfpmean=mean(bfp, na.rm=TRUE), 
#                                                                    timemean=mean(time))

data = data %>% group_by(root, clade, div, episode) %>% mutate(timemean=mean(time))
datas = data %>% group_by(root, clade, div, episode, timemean) %>% summarise(ppmax=max(pp), 
                                                                             ppmean=mean(pp),
                                                                             bfpmax=max(bf), 
                                                                             bfmean=mean(bf))

#---- Beautifying tables for plotting --------------------------------------------------------------

geos = subset(geo, time < max(data$timemean))
geoss = subset(era, time < max(data$timemean))

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

data$colour = factor(data$colour, levels=c("royalblue1", "steelblue2", "steelblue4",
                                           "forestgreen", "orange2",
                                           "yellow1", "hotpink1", "darkseagreen3",
                                           "darkorchid2", "darkorchid3" , "darkorchid4"
))

data$root = factor(data$root, levels=c("rootD", "rootA"))
data$clade = factor(data$clade, levels=c("Amoebozoa", "Nucletmycea", "Holozoa",
                                         "Metamonada", "Discoba",
                                         "Haptista", "Cryptista", "Archaeplastida",
                                         "Rhizaria", "Stramenopila", "Alveolata"
))
datas$clade = factor(datas$clade, levels=c("Amoebozoa", "Nucletmycea", "Holozoa",
                                           "Metamonada", "Discoba",
                                           "Haptista", "Cryptista", "Archaeplastida",
                                           "Rhizaria", "Stramenopila", "Alveolata"
))

#---- Plotting -------------------------------------------------------------------------------------

pdf(paste0(files$outPrefix, "_individualClades.pdf"), width=11.69, height=8.27, paper='special')
# pdf(paste0(files$outPrefix, "_individualClades_BF.pdf"), width=11.69, height=8.27, paper='special')
for(cladei in sort(unique(datas$clade))){
  cat("\r  Plotting clade ", cladei, "\t", grep(cladei, sort(unique(datas$clade))), "/", length(unique(datas$clade)), "          ", sep="")
  datass = subset(data, clade==cladei)
  datasss = subset(datas, clade==cladei)
  geos = subset(geo, time < max(datass$timemean))
  geoss = subset(era, time < max(datass$timemean))
  plotME = ggplot(datass) + 
    geom_vline(xintercept = geos$time, color="grey80") +
    geom_vline(xintercept = geoss$time, color="grey60") +
    geom_bar(aes(x = timemean, y = pp, fill=clade), alpha=0.2, stat="identity", position=position_dodge()) +
    geom_point(data=datasss, aes(x = timemean, y = ppmean, colour=clade), alpha=1, colour="black") +
    # geom_bar(aes(x = timemean, y = bf, fill=clade), alpha=0.2, stat="identity", position=position_dodge()) +
    # geom_point(data=datasss, aes(x = timemean, y = bfmean, colour=clade), alpha=1, colour="black") +
    facet_grid(root~div, scales="free")+
    scale_x_reverse(breaks=seq(0, max(datass$timemean), 100)) + 
    labs(title=paste("Mass extinctions", cladei),
         x="Million years ago",
         y="Posterior probabilities of Mass Extinctions") +
         # y="Bayes Factor for Mass Extinctions") +
    theme_classic() + 
    scale_fill_manual(values=as.character(sort(unique(datass$colour)))) +
    theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA))
  plot(plotME)
}; rm(cladei, datass, geos, geoss); cat("\rDone\n")
dev.off()

#----
#---- Summarizing tables per period and clade ------------------------------------------------------

data =  data %>% group_by(root, clade, episode) %>% mutate(timemean=mean(time))

data$period = NA
for(p in geo$period){
  t = subset(geo, period==p)$time
  data$period[which(data$timemean < t)] = p
};rm(p, t)
any(is.na((data$period)))
data$period = ifelse(data$period=="Paleogene" | data$period=="Neogene" | data$period=="Quaternary", "Cenozoic", data$period)

datap = data %>% group_by(root, clade, div, period) %>% summarise(ppmax=max(pp), ppmean=mean(pp),
                                                                  bfmax=max(bf), bfmean=mean(bf))

#---- Beautifying tables for plotting --------------------------------------------------------------

datap$clade = factor(datap$clade, levels=c("Amoebozoa", "Nucletmycea", "Holozoa",
                                           "Metamonada", "Discoba",
                                           "Haptista", "Cryptista", "Archaeplastida",
                                           "Rhizaria", "Stramenopila", "Alveolata"))

{datap$clades = as.character(datap$clade)
  datap$clades[which(datap$clade=="Amoebozoa")] = "Am"
  datap$clades[which(datap$clade=="Nucletmycea")] = "N"
  datap$clades[which(datap$clade=="Holozoa")] = "Ho"
  datap$clades[which(datap$clade=="Metamonada")] = "M"
  datap$clades[which(datap$clade=="Discoba")] = "D"
  datap$clades[which(datap$clade=="Haptista")] = "Ha"
  datap$clades[which(datap$clade=="Cryptista")] = "C"
  datap$clades[which(datap$clade=="Archaeplastida")] = "Ar"
  datap$clades[which(datap$clade=="Rhizaria")] = "R"
  datap$clades[which(datap$clade=="Stramenopila")] = "S"
  datap$clades[which(datap$clade=="Alveolata")] = "Al"
  datap$clades = factor(datap$clades, levels=c("Am","N","Ho","M","D","Ha","C","Ar","R","S","Al"))}

{datap$colour = as.character(datap$clade)
  datap$colour[which(datap$colour=="Amoebozoa")]="royalblue1"
  datap$colour[which(datap$colour=="Nucletmycea")]="steelblue2"
  datap$colour[which(datap$colour=="Holozoa")]="steelblue4"
  datap$colour[which(datap$colour=="Metamonada")]="forestgreen"
  datap$colour[which(datap$colour=="Discoba")]="orange2"
  datap$colour[which(datap$colour=="Haptista")]="yellow1"
  datap$colour[which(datap$colour=="Cryptista")]="hotpink1"
  datap$colour[which(datap$colour=="Archaeplastida")]="darkseagreen3"
  datap$colour[which(datap$colour=="Rhizaria")]="darkorchid2"
  datap$colour[which(datap$colour=="Stramenopila")]="darkorchid3"
  datap$colour[which(datap$colour=="Alveolata")]="darkorchid4"
  datap$colour = factor(datap$colour,
                        levels=c("royalblue1","steelblue2","steelblue4",
                                 "forestgreen","orange2",
                                 "yellow1","hotpink1","darkseagreen3",
                                 "darkorchid2", "darkorchid3","darkorchid4"))}

datap$root = factor(datap$root, levels=c("rootD", "rootA"))


datap$period = factor(datap$period, levels=c("Statherian",
                                             "Calymnian", "Ectasian", "Stenian", 
                                             "Tonian", "Cryogenian", "Ediacaran", 
                                             "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", 
                                             "Triassic", "Jurassic", "Cretaceous", 
                                             "Cenozoic"))

#---- Plotting -------------------------------------------------------------------------------------

(plotMEsum = ggplot(datap)+
   # geom_bar(aes(x=clades, y=ppmax, fill=clade), alpha=0.2, stat="identity", position=position_dodge())+
   # geom_bar(aes(x=clades, y=ppmean, fill=clade), stat="identity", position=position_dodge())+
   geom_bar(aes(x=clades, y=bfmax, fill=clade), alpha=0.2, stat="identity", position=position_dodge())+
   geom_bar(aes(x=clades, y=bfmean, fill=clade), stat="identity", position=position_dodge())+
   facet_grid(root+div~period, scales="free")+
   scale_fill_manual(values=as.character(sort(unique(datap$colour)))) +
   theme_minimal()+
   # labs(x=NULL, y="Summary of Posterior Probabilities")+
   labs(x=NULL, y="Summary of 2 Ln(Bayes Factor)")+
   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
         legend.position="bottom", legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=6), legend.title = element_blank())+ 
   guides(fill = guide_legend(nrow = 1)))

pdf(paste0(files$outPrefix, "_summaryPeriods_BF.pdf"), width=11.69, height=8.27, paper='special'); plot(plotMEsum); dev.off()

#----
#----
