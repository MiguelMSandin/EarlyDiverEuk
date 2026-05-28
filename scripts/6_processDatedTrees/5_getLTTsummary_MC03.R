#----
#---- loading packages and data.frames -------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----
#---- set working directory and file names ---------------------------------------------------------
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/")

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirsMC1="stepDating/MC01-2/", # The directory where the MC01 dated trees are,
             dirsMC2="stepDating/MC02-2/", # The directory where the MC02 dated trees are
             dirsMC3="stepDating/MC03-2/", # The directory where the MC03 dated trees are
             dirsMC4="stepDating/MC04/", # The directory where the MC04 dated trees are
             extension="LTT-data.tsv", # The common extension between all the LTT tables
             factor=-1, # A factor to multiply the branch lengths for plotting in millions of years ago
             out="stepDating/MC_summary_cleaned_LTT") # The output files without extension for the tsv table (*.tsv) and the pdf (*.pdf)

#----
#---- Read tables ----------------------------------------------------------------------------------

# First search all the tables in the directories recursively
(files$tablesMC1 = grep("_LTT-data.tsv", dir(files$dirsMC1, recursive=TRUE), value=TRUE))
(files$tablesMC2 = grep("_LTT-data.tsv", dir(files$dirsMC2, recursive=TRUE), value=TRUE))
(files$tablesMC3 = grep("_LTT-data.tsv", dir(files$dirsMC3, recursive=TRUE), value=TRUE))
(files$tablesMC4 = grep("_LTT-data.tsv", dir(files$dirsMC4, recursive=TRUE), value=TRUE))

# Now concatenate all tables into one
data = data.frame()
i = 0
for(file in files$tablesMC1){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC1), sep="")
	f = fread(paste0(files$dirsMC1, "/", file))
	f = subset(f, tree=="Main_tree" | tree=="Stramenopila" | tree=="Alveolata" | tree=="Archaeplastida")
	f$calibration = "MC01"
	f$file = file %>% gsub("root.\\/", "", .) %>% gsub("/.*", "", .)
	f$root = gsub("/.*", "", file)
	data = rbind(data, f)
}; cat("\n")

i = 0
for(file in files$tablesMC2){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC2), sep="")
	f = fread(paste0(files$dirsMC2, "/", file))
	f = subset(f, tree=="Main_tree" | tree=="Archaeplastida" | tree=="Stramenopila" | tree=="Alveolata")
	f$calibration = "MC02"
	f$file = file %>% gsub("root.\\/", "", .) %>% gsub("/.*", "", .)
	f$root = gsub("/.*", "", file)
	data = rbind(data, f)
}; rm(i, file, f); cat("\n")

i = 0
for(file in files$tablesMC3){
  cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC3), sep="")
  f = fread(paste0(files$dirsMC3, "/", file))
  f = subset(f, tree=="Main_tree" | tree=="Stramenopila")
  f$calibration = "MC03"
  f$file = file %>% gsub("root.\\/", "", .) %>% gsub("/.*", "", .)
  f$root = gsub("/.*", "", file)
  data = rbind(data, f)
}; rm(i, file, f); cat("\n")


i = 0
for(file in files$tablesMC4){
  cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC4), sep="")
  f = fread(paste0(files$dirsMC4, "/", file))
  f = subset(f, tree=="Main_tree" | tree=="Archaeplastida")
  f$calibration = "MC04"
  f$file = file %>% gsub("root.\\/", "", .) %>% gsub("/.*", "", .)
  f$root = gsub("/.*", "", file)
  data = rbind(data, f)
}; rm(i, file, f); cat("\n")


data

unique(data$tree)
unique(data$calibration)
unique(data$root)
unique(data$file)

# ltt.plot function counts nodes and not lineages, so we have to transform nodes (n) to lineages (n+1)
data$N = data$N + 1
data$logN = log(data$N)

#----
#---- Get stats ------------------------------------------------------------------------------------

# ancovaCal = aov(time ~ logN * calibration, data)
# ancovaFile = aov(time ~ logN * file, data)
# 
# summary(ancovaCal)
# summary(ancovaFile)

subset(data, N==2) %>% group_by(root, tree, calibration) %>% 
  summarise(mean=mean(time), sd=sd(time), median=median(time), min=min(time), max=max(time))

tmp = subset(data, N==2)
summary(lm(time~tree*root*calibration, tmp))

tmp = subset(data, tree=="Archaeplastida" & N==2)
ggplot(tmp, aes(x=time, y=paste0(calibration, "-", root)))+geom_boxplot()+theme_bw()+labs(title=paste(unique(tmp$tree), "root"))
summary(lm(time~root*calibration, tmp))

tmp = subset(data, tree=="Stramenopila" & N==2)
ggplot(tmp, aes(x=time, y=paste0(calibration, "-", root)))+geom_boxplot()+theme_bw()+labs(title=paste(unique(tmp$tree), "root"))
summary(lm(time~root*calibration, tmp))

tmp = subset(data, tree=="Main_tree" & N==2)
ggplot(tmp, aes(x=time, y=paste0(calibration, "-", root)))+geom_boxplot()+theme_bw()+labs(title=paste(unique(tmp$tree), "root"))
summary(lm(time~root*calibration, tmp))




datas = subset(data, N==2) %>% group_by(root, tree, calibration) %>% 
  summarise(mean=mean(time), sd=sd(time), median=median(time), min=min(time), max=max(time))



# Summarise table
# data$time = data$time * -1
datas = data %>% group_by(calibration, tree, root, N) %>% summarise(median=quantile(time, 0.50),
                                                                    q05=quantile(time, 0.05),
                                                                    q95=quantile(time, 0.95),
                                                                    replicates=length(unique(file)))

# only use lineages and time points with at least 50% of the replicates
datas = subset(datas, (calibration=="MC01" & replicates >= 16) | ((calibration=="MC02" | calibration=="MC03" | calibration=="MC04") & replicates >= 5))

#----
#---- beautifying and plotting the summary dataset for plotting ------------------------------------

sort(unique(datas$tree))
sort(unique(datas$calibration))

datas$treeCal = paste0(datas$tree, "-", datas$calibration)
sort(unique(datas$treeCal))
datas$treeCal = factor(datas$treeCal,
                       levels=c("Main_tree-MC01", "Main_tree-MC02", "Main_tree-MC03", "Main_tree-MC04",
                                "Archaeplastida-MC01", "Archaeplastida-MC02", "Archaeplastida-MC04",
                                "Stramenopila-MC01", "Stramenopila-MC02", "Stramenopila-MC03"))

{datas$colour = as.character(datas$treeCal)
  datas$colour[which(datas$colour=="Main_tree-MC01")]="black"
  datas$colour[which(datas$colour=="Main_tree-MC02")]="grey20"
  datas$colour[which(datas$colour=="Main_tree-MC03")]="grey40"
  datas$colour[which(datas$colour=="Main_tree-MC04")]="grey60"
  datas$colour[which(datas$colour=="Archaeplastida-MC01")]="darkseagreen4"
  datas$colour[which(datas$colour=="Archaeplastida-MC02")]="darkseagreen3"
  datas$colour[which(datas$colour=="Archaeplastida-MC04")]="darkseagreen1"
  datas$colour[which(datas$colour=="Stramenopila-MC01")]="#680D96"
  datas$colour[which(datas$colour=="Stramenopila-MC02")]="#8611C0"
  datas$colour[which(datas$colour=="Stramenopila-MC03")]="#A314EB"}
datas$colour = factor(datas$colour, levels=c("black", "grey20", "grey40", "grey60", "grey80", 
                                             "#B23AEE", "#C469F2", "#D593F6", 
                                             "darkseagreen4", "darkseagreen3", "darkseagreen1", 
                                             "#680D96", "#8611C0", "#A314EB"))

datas$tree = factor(datas$tree, levels=c("Main_tree", "Archaeplastida", "Alveolata", "Stramenopila"))

datas$root = factor(datas$root, levels=c("rootD", "rootA"))

geos = subset(geo, time > min(datas$q05, na.rm=TRUE))

(lttplot = ggplot(datas)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_ribbon(aes(xmin=q95, xmax=q05, y=N, fill=colour, group=treeCal), alpha=0.2, colour = NA) +
    geom_line(aes(x=median, y=N, colour=colour, group=treeCal)) +
    # geom_text(data=geos, aes(x=mid, y=max(datas$N)+5, label=period), size=2)+
    facet_wrap(~root, nrow=2, ncol=1)+
    scale_x_continuous(breaks=seq((round(min(datas$q05, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(datas$q05, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_fill_manual(values=as.character(sort(unique(datas$colour))),
                      labels=as.character(sort(unique(datas$treeCal))))+
    scale_color_manual(values=as.character(sort(unique(datas$colour))),
                       labels=as.character(sort(unique(datas$treeCal))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))

pdf(paste0(files$out, "_summary.pdf"), width=11.69, height=8.27*2, paper='special'); plot(lttplot); dev.off()

#----
