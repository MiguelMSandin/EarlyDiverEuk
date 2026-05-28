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
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dir="./", # The directory where the MC01 dated trees are,
             extension="LTT-data.tsv", # The common extension between all the LTT tables
             factor=-1, # A factor to multiply the branch lengths for plotting in millions of years ago
             out="plots/LTT/LTT_MC01_all_R1") # The output files without extension for the tsv table (*.tsv) and the pdf (*.pdf)

if(!dir.exists(dirname(files$out))){dir.create(dirname(files$out))}

#----
#---- Read tables ----------------------------------------------------------------------------------

# First search all the tables in the directories recursively
(files$tables = grep("_LTT-data.tsv", dir(files$dir, recursive=TRUE), value=TRUE))

# Now concatenate all tables into one
data = data.frame()
i = 0
for(file in files$tables){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tables), sep="")
	f = fread(paste0(files$dir, "/", file))
	f$root = gsub("\\/.*", "", file)
	f$file = file %>% gsub("root.\\/", "", .) %>% gsub("\\/step3.*", "", .)
	data = rbind(data, f)
};rm(i, f, file); cat("\n")

# ltt.plot function counts nodes and not lineages, so we have to transform nodes (n) to lineages (n+1)
data$N = data$N + 1
data$logN = log(data$N)

#----
#---- Get stats ------------------------------------------------------------------------------------

# test for differences in the rooting scenarios
tmp = subset(data, N == 2 & tree =="Main_tree")

ggplot(tmp) + geom_density(aes(x=time, fill=root, alpha=0.2))
var.test(time ~ root, data=tmp)
# t.test(median ~ root, data=tmp2)
(aov(time ~ root, data=tmp))
anova((aov(time ~ root, data=tmp)))

# Now check for differences within the supergroups
tmp = unique(subset(data, N >= 300)$tree)
tmp = subset(data, N == 2 & tree %in% tmp)
pdf(paste0(files$out, "_differences.pdf"), width=11.69, height=8.27, paper='special')
for(group in unique(tmp$tree)){
	ss = subset(tmp, tree == group)
	v = var.test(time ~ root, data=ss)
	if(v$p.value > 0.05){
		a = aov(time ~ root, data=ss)
		a = anova(a)
		test = "ANOVA"
		p = a$`Pr(>F)`
	}else{
		t = t.test(time ~ root, data=ss)
		test = "T-test"
		p = t$p.value
	}
	s = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", ifelse(p < 1, " ", )))))
	plot = ggplot(ss) + 
		geom_density(aes(x=time, fill=root, alpha=0.2)) + 
		labs(title=group, subtitle=paste0(test, ": ", round(p, 4), " ", s))+
		theme_minimal()
	plot(plot)
}; rm(group, ss, v, a, t, test, p, s); dev.off()


# Check average appearance times
tmp = subset(data, N == 2)

tmp2 = tmp %>% group_by(root, tree) %>% summarise(median=median(time), 
												  min=min(time),
												  max=max(time))
(tmp2 = tmp2[order(tmp2$median),])
tmp$tree = factor(tmp$tree, levels=unique(tmp2$tree))

ggplot(tmp) +
	geom_point(aes(y=tree, x=time, colour=root))
ggplot(tmp, aes(y=tree, x=time, colour=root)) +
	geom_boxplot()+
	geom_jitter(alpha=0.2)+
	scale_x_continuous(breaks=seq(round(min(tmp$time), -3), 0, 100), limits=c(round(min(tmp$time), -3), 0)) +
	theme_bw()

#----
#---- Summarise, subset and export table -----------------------------------------------------------

# Summarise table
data$time = data$time * files$factor
summ = data %>% group_by(root, tree, N, logN) %>% summarise(median=quantile(time, 0.50),
                                                            # min=min(time),
                                                            q05=quantile(time, 0.05),
                                                            q95=quantile(time, 0.95),
                                                            # max=max(time),
                                                            replicates=length(unique(file)))

# only use lineages and time points with at least 50% of the replicates
tmp = max(summ[!grepl("MC01", summ$tree),]$replicates)
summ = subset(summ, replicates >= tmp*0.5)

# Only use groups with at least 300 lineages
summ = summ[summ$tree %in% names(table(summ$tree)[table(summ$tree)/2 > 300]),]

# Exporting the dataset ____________________________________________________________________________
write.table(summ, paste0(files$out, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
#---- Reading table if already available -----------------------------------------------------------

summ = fread(paste0(files$out, ".tsv"))

#----
#---- beautifying the dataset for plotting ---------------------------------------------------------

sort(unique(summ$root))
summ$root = factor(summ$root, levels=c("rootD", "rootA"))

sort(unique(summ$tree))
summ$tree = factor(summ$tree, 
                   levels=c("Main_tree",
                            "Amoebozoa", "Nucletmycea", "Holozoa", 
                            "Metamonada", "Discoba", 
                            "Haptista", "Cryptista", "Archaeplastida", 
                            "Rhizaria", "Stramenopila", "Alveolata"))

{summ$colour = as.character(summ$tree)
	summ$colour[which(summ$colour=="Main_tree")]="black"
	summ$colour[which(summ$colour=="Main_tree")]="grey80"
    summ$colour[which(summ$colour=="Amoebozoa")]="royalblue1"
    summ$colour[which(summ$colour=="Apusomonadida")]="lightskyblue1"
    # summ$colour[which(summ$colour=="Breviatea")]="lightskyblue2"
    summ$colour[which(summ$colour=="Nucletmycea")]="steelblue2"
    summ$colour[which(summ$colour=="Holozoa")]="steelblue4"
    # summ$colour[which(summ$colour=="Ancyromonadida")]="grey20"
    # summ$colour[which(summ$colour=="CRuMs")]="aquamarine1"
    summ$colour[which(summ$colour=="Metamonada")]="forestgreen"
    summ$colour[which(summ$colour=="Discoba")]="orange2"
    summ$colour[which(summ$colour=="Haptista")]="yellow1"
    summ$colour[which(summ$colour=="Cryptista")]="hotpink1"
    # summ$colour[which(summ$colour=="Hemimastigophora")]="grey80"
    summ$colour[which(summ$colour=="Archaeplastida")]="darkseagreen3"
    # summ$colour[which(summ$colour=="Telonemia")]="plum4"
    summ$colour[which(summ$colour=="Rhizaria")]="darkorchid2"
    summ$colour[which(summ$colour=="Stramenopila")]="darkorchid3"
    summ$colour[which(summ$colour=="Alveolata")]="darkorchid4"}

summ$colour = factor(summ$colour, 
                      levels=c("black", "grey80",
                               "royalblue1", "lightskyblue1", "steelblue2", "steelblue4",
                               "forestgreen", "orange2",
                               "yellow1", "hotpink1", "darkseagreen3",
                               "darkorchid2", "darkorchid3", "darkorchid4"))
sort(unique(summ$colour))

#----
#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table
geos = subset(geo, time > min(-summ$q95, na.rm=TRUE))

(lttplot = ggplot(summ)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    geom_line(aes(x=-median, y=N, color=tree)) +
    facet_wrap(~root, ncol=2)+
    # geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Lineages", x="Time (Ma)")+
    theme(legend.position = "none"))

pdf(paste0(files$out, ".pdf"), width=11.69*1.5, height=8.27*0.8, paper='special'); plot(lttplot); dev.off()

(lttplotD = ggplot(subset(summ, root=="rootD"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    geom_line(aes(x=-median, y=N, color=tree)) +
    geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Lineages", x="Time (Ma)"))
pdf(paste0(files$out, "_rootD.pdf"), width=11.69, height=8.27, paper='special'); plot(lttplotD); dev.off()

(lttplotA = ggplot(subset(summ, root=="rootA"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    geom_line(aes(x=-median, y=N, color=tree)) +
    geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Lineages", x="Time (Ma)"))
pdf(paste0(files$out, "_rootA.pdf"), width=11.69, height=8.27, paper='special'); plot(lttplotA); dev.off()


#----
