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
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/")

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirsMC1="stepDating/MC01-2/", # The directory where the MC01 dated trees are,
             dirsMC2="stepDating/MC02-2/", # The directory where the MC02 dated trees are
             dirsMC3="stepDating/MC03-2/", # The directory where the MC03 dated trees are
             extension="LTT-data.tsv", # The common extension between all the LTT tables
             factor=-1, # A factor to multiply the branch lengths for plotting in millions of years ago
             out="stepDating/MC03-2/plots/MC03_cleaned_LTT") # The output files without extension for the tsv table (*.tsv) and the pdf (*.pdf)

#----
#---- Read tables ----------------------------------------------------------------------------------

# First search all the tables in the directories recursively
files$tablesMC1 = grep("_LTT-data.tsv", dir(files$dirsMC1, recursive=TRUE), value=TRUE)
files$tablesMC2 = grep("_LTT-data.tsv", dir(files$dirsMC2, recursive=TRUE), value=TRUE)
files$tablesMC3 = grep("_LTT-data.tsv", dir(files$dirsMC3, recursive=TRUE), value=TRUE)

files

# Now concatenate all tables into one
data = data.frame()
i = 0
for(file in files$tablesMC1){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC1), sep="")
	f = fread(paste0(files$dirsMC1, "/", file))
	f = subset(f, tree=="Main_tree" | tree=="Stramenopila")
	f$calibration = "MC01"
	f$file = file %>% gsub("root.\\/", "", .) %>% gsub("/.*", "", .)
	f$root = gsub("/.*", "", file)
	data = rbind(data, f)
}; cat("\n")

i = 0
for(file in files$tablesMC2){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC2), sep="")
	f = fread(paste0(files$dirsMC2, "/", file))
	f = subset(f, tree=="Main_tree" | tree=="Stramenopila")
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

data

unique(data$tree)
unique(data$calibration)
unique(data$root)
unique(data$file)

#----
#---- Get stats ------------------------------------------------------------------------------------

# ancovaCal = aov(time ~ logN * calibration, data)
# ancovaFile = aov(time ~ logN * file, data)
# 
# summary(ancovaCal)
# summary(ancovaFile)

subset(data, N==1) %>% group_by(root, tree, calibration) %>% 
  summarise(mean=mean(time), sd=sd(time), median=median(time), min=min(time), max=max(time))

tmp = subset(data, N==1)
summary(lm(time~tree*root*calibration, tmp))

tmp = subset(data, tree=="Stramenopila" & N==1)
summary(lm(time~root*calibration, tmp))

tmp = subset(data, tree=="Main_tree" & N==1)
summary(lm(time~root*calibration, tmp))




datas = data %>% group_by(root, tree, calibration) %>% 
  summarise(mean=mean(time), sd=sd(time), median=median(time), min=min(time), max=max(time))



# Summarise table
# data$time = data$time * -1
datas = data %>% group_by(calibration, tree, root, N) %>% summarise(median=quantile(time, 0.50),
                                                                    q05=quantile(time, 0.05),
                                                                    q95=quantile(time, 0.95),
                                                                    replicates=length(unique(file)))

# only use lineages and time points with at least 50% of the replicates
datas = subset(datas, (calibration=="MC01" & replicates >= 16) | ((calibration=="MC02" | calibration=="MC03") & replicates >= 5))

#----
#---- beautifying the dataset for plotting ---------------------------------------------------------

sort(unique(data$tree))
sort(unique(data$calibration))

data$treeCal = paste0(data$tree, "-", data$calibration)

data$treeCal = factor(data$treeCal, 
                      levels=c("Main_tree-MC01", "Main_tree-MC02", "Main_tree-MC03", 
                               "Stramenopila-MC01", "Stramenopila-MC02", "Stramenopila-MC03"))

data$root = factor(data$root, levels=c("rootD", "rootA"))

{data$colour = as.character(data$treeCal)
	data$colour[which(data$colour=="Main_tree-MC01")]="grey60"
	data$colour[which(data$colour=="Main_tree-MC02")]="grey90"
	data$colour[which(data$colour=="Main_tree-MC03")]="black"
	data$colour[which(data$colour=="Stramenopila-MC01")]="#D5ABEA"
	data$colour[which(data$colour=="Stramenopila-MC02")]="#E4C9F1"
	data$colour[which(data$colour=="Stramenopila-MC03")]="darkorchid3"}

data$colour = factor(data$colour, 
					 levels=c("grey60", "grey90", "black", 
					 		 "#D5ABEA", "#E4C9F1", "darkorchid3"))

data$group = paste0(data$file, "-", data$calibration, "-" , data$tree)

#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table
# datas = subset(data, file=="fRAcat2-2" | file=="fRAgamma1-2")
# datas = subset(data, tree=="Stramenopila")

geos = subset(geo, time > min(data$time, na.rm=TRUE))

(lttplot = ggplot(data)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_line(aes(x=time, y=N, group=group, colour=treeCal)) +
    geom_text(data=geos, aes(x=mid, y=max(data$N)+5, label=period), size=2)+
    facet_wrap(~root, nrow=2, ncol=1)+
    scale_x_continuous(breaks=seq((round(min(data$time, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(data$time, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_color_manual(values=as.character(sort(unique(data$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))

pdf(paste0(files$out, ".pdf"), width=11.69, height=8.27*2, paper='special'); plot(lttplot); dev.off()

#----
#---- beautifying and plotting the summary dataset for plotting ------------------------------------

sort(unique(datas$tree))
sort(unique(datas$calibration))

datas$treeCal = paste0(datas$tree, "-", datas$calibration)

datas$treeCal = factor(datas$treeCal, 
                      levels=c("Main_tree-MC01", "Main_tree-MC02", "Main_tree-MC03", 
                               "Stramenopila-MC01", "Stramenopila-MC02", "Stramenopila-MC03"))

datas$root = factor(datas$root, levels=c("rootD", "rootA"))

{datas$colour = as.character(datas$treeCal)
  datas$colour[which(datas$colour=="Main_tree-MC01")]="grey60"
  datas$colour[which(datas$colour=="Main_tree-MC02")]="grey90"
  datas$colour[which(datas$colour=="Main_tree-MC03")]="black"
  datas$colour[which(datas$colour=="Stramenopila-MC01")]="#D5ABEA"
  datas$colour[which(datas$colour=="Stramenopila-MC02")]="#E4C9F1"
  datas$colour[which(datas$colour=="Stramenopila-MC03")]="darkorchid3"}

datas$colour = factor(datas$colour, 
                     levels=c("grey60", "grey90", "black", 
                              "#D5ABEA", "#E4C9F1", "darkorchid3"))

geos = subset(geo, time > min(datas$q05, na.rm=TRUE))

(lttplot = ggplot(datas)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_ribbon(aes(xmin=q95, xmax=q05, y=N, fill=treeCal), alpha=0.2, colour = NA) +
    geom_line(aes(x=median, y=N, group=group, colour=treeCal)) +
    geom_text(data=geos, aes(x=mid, y=max(data$N)+5, label=period), size=2)+
    facet_wrap(~root, nrow=2, ncol=1)+
    scale_x_continuous(breaks=seq((round(min(data$time, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(data$time, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_fill_manual(values=as.character(sort(unique(data$colour))))+
    scale_color_manual(values=as.character(sort(unique(data$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))

pdf(paste0(files$out, "_summary.pdf"), width=11.69, height=8.27*2, paper='special'); plot(lttplot); dev.off()

#----
