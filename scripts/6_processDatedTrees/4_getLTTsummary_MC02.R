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
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/")

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(dirsMC1="MC01-2/", # The directory where the MC01 dated trees are,
             dirsMC2="MC02-2/", # The directory where the MC02 dated trees are
             extension="LTT-data.tsv", # The common extension between all the LTT tables
             factor=-1, # A factor to multiply the branch lengths for plotting in millions of years ago
             out="MC02-2/plots/LTT_MC02") # The output files without extension for the tsv table (*.tsv) and the pdf (*.pdf)

if(!dir.exists(dirname(files$out))){dir.create(dirname(files$out))}

#----
#---- Read tables ----------------------------------------------------------------------------------

# First search all the tables in the directories recursively
files$tablesMC1 = grep("_LTT-data.tsv", dir(files$dirsMC1, recursive=TRUE), value=TRUE)
files$tablesMC2 = grep("_LTT-data.tsv", dir(files$dirsMC2, recursive=TRUE), value=TRUE)

files

# Now concatenate all tables into one
data = data.frame()
i = 0
for(file in files$tablesMC1){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC1), sep="")
	f = fread(paste0(files$dirsMC1, "/", file))
	f = subset(f, tree=="Main_tree")
	f$tree = "MC01-main_tree"
	f$root = gsub("\\/.*", "", file)
	f$file = gsub(".*\\/", "MC01-", file)
	data = rbind(data, f)
}; cat("\n")

i = 0
for(file in files$tablesMC2){
	cat("\r  Reading files ", (i = i + 1), "/", length(files$tablesMC2), sep="")
	f = fread(paste0(files$dirsMC2, "/", file))
	f$root = gsub("\\/.*", "", file)
	f$file = gsub(".*\\/", "MC02-", file)
	data = rbind(data, f)
}; rm(i, file, f); cat("\n")

#----
#---- Get stats ------------------------------------------------------------------------------------

# Get some sort of comparisons for the early diversification

tmp = subset(data, (tree=="MC01-main_tree" | tree=="Main_tree") & N <= 100)
tmp = tmp %>% group_by(tree, root, N) %>% summarise(median=median(time), 
                                                    min=min(time),
                                                    max=max(time))
ggplot(tmp, aes(y=N))+
	geom_point(aes(x=median))+
	geom_point(aes(x=min), colour="grey80")+
	geom_point(aes(x=max), colour="grey80")+theme_minimal()

lm(logN ~ time, subset(data, tree=="MC01-main_tree" & N <= 20))$coeff[2]
lm(logN ~ time, subset(data, tree=="Main_tree" & N <= 20))$coeff[2]

tmp = subset(data, (tree=="MC01-main_tree" | tree=="Main_tree") & (N == 1 | N == 2 | N == 4 | N == 6 | N == 8 | N == 10 |  N == 20 | N == 50 | N == 100))
tmp = tmp %>% group_by(tree, root, N) %>% summarise(median=median(time), min=min(time), max=max(time))
for(l in c(2, 4, 6, 8, 10, 20, 50, 100)){
	cat("MC01 reaches ", l, " lineages in ",
		-round(subset(tmp, tree=="MC01-main_tree" & N == 1)$median - subset(tmp, tree=="MC01-main_tree" & N == l)$median), " (", 
		-round(subset(tmp, tree=="MC01-main_tree" & N == 1)$min - subset(tmp, tree=="MC01-main_tree" & N == l)$min), "-", 
		-round(subset(tmp, tree=="MC01-main_tree" & N == 1)$max - subset(tmp, tree=="MC01-main_tree" & N == l)$max), 
		") years, and MC02 in ",
		-round(subset(tmp, tree=="Main_tree" & N == 1)$median - subset(tmp, tree=="Main_tree" & N == l)$median), " (", 
		-round(subset(tmp, tree=="Main_tree" & N == 1)$min - subset(tmp, tree=="Main_tree" & N == l)$min), "-", 
		-round(subset(tmp, tree=="Main_tree" & N == 1)$max - subset(tmp, tree=="Main_tree" & N == l)$max), 
		") years, about ",
		-round((subset(tmp, tree=="MC01-main_tree" & N == 1)$median - subset(tmp, tree=="MC01-main_tree" & N == l)$median)-
			   	(subset(tmp, tree=="Main_tree" & N == 1)$median - subset(tmp, tree=="Main_tree" & N == l)$median)),
		" years of difference.",
		sep="", end="\n")
}; rm(l)


# Plot origin times
tmp = subset(data, (tree=="MC01-main_tree" | tree=="Main_tree") & N == 1)
tmp %>% group_by(tree, root) %>% summarise(median=median(time), 
                                           min=min(time),
                                           max=max(time))
ggplot(tmp) +
	geom_point(aes(y=tree, x=time))

# Get slopes per 100 years of every supergroup and replicate
slopes = data.frame()
for(group in unique(data$tree)){
    i = 0
    tmp = data[data$tree == group,]
    if(max(tmp$N) >= 300){
        for(replicate in unique(tmp$file)){
            cat("\r  Analysing ", group, " (", (i = i + 1), "/", length(unique(tmp$file)), ")          ", sep="")
            tmp1 = tmp[tmp$file == replicate,]
            M = round(min(tmp1$time))+100
            j = 0
            for(m in c(seq(M, 0, 100), 0)){
                j = j + 100
                tmp2 = subset(tmp1, time <= m)
                slopes = rbind(slopes, data.frame(group=group, file=replicate, 
                                                   age=m, meassured = j,
                                                   N=max(tmp2$N), logN=max(tmp2$logN),
                                                   slope=lm(tmp2$logN ~ tmp2$time)$coeff[2]))
            }
        }
    }else{
        cat("  ", group, " not analysed: less than 300 lineages in tree.", sep="")
    }
    cat("\n")
}; rm(i, group, replicate, tmp, tmp1, M, m, tmp2, j); cat("Done\n")

sort(unique(slopes$group))
slopes$group = factor(slopes$group, 
                    levels=c("Main_tree", "MC01-main_tree", 
                             "Amoebozoa", "Nucletmycea", "Holozoa", 
                             "Metamonada", "Discoba", 
                             "Haptista", "Cryptista", "Archaeplastida", 
                             "Rhizaria", "Stramenopila", "Alveolata"))

{slopes$colour = as.character(slopes$group)
	slopes$colour[which(slopes$colour=="Main_tree")]="black"
	slopes$colour[which(slopes$colour=="MC01-main_tree")]="grey80"
    slopes$colour[which(slopes$colour=="Amoebozoa")]="royalblue1"
    slopes$colour[which(slopes$colour=="Apusomonadida")]="lightskyblue1"
    # slopes$colour[which(slopes$colour=="Breviatea")]="lightskyblue2"
    slopes$colour[which(slopes$colour=="Nucletmycea")]="steelblue2"
    slopes$colour[which(slopes$colour=="Holozoa")]="steelblue4"
    # slopes$colour[which(slopes$colour=="Ancyromonadida")]="grey20"
    # slopes$colour[which(slopes$colour=="CRuMs")]="aquamarine1"
    slopes$colour[which(slopes$colour=="Metamonada")]="forestgreen"
    slopes$colour[which(slopes$colour=="Discoba")]="orange2"
    slopes$colour[which(slopes$colour=="Haptista")]="yellow1"
    slopes$colour[which(slopes$colour=="Cryptista")]="hotpink1"
    # slopes$colour[which(slopes$colour=="Hemimastigophora")]="grey80"
    slopes$colour[which(slopes$colour=="Archaeplastida")]="darkseagreen3"
    # slopes$colour[which(slopes$colour=="Telonemia")]="plum4"
    slopes$colour[which(slopes$colour=="Rhizaria")]="darkorchid2"
    slopes$colour[which(slopes$colour=="Stramenopila")]="darkorchid3"
    slopes$colour[which(slopes$colour=="Alveolata")]="darkorchid4"}

slopes$colour = factor(slopes$colour, 
                      levels=c("black", "grey80",
                               "royalblue1", "lightskyblue1", "steelblue2", "steelblue4",
                               "forestgreen", "orange2",
                               "yellow1", "hotpink1", "darkseagreen3",
                               "darkorchid2", "darkorchid3", "darkorchid4"))

slopes$meassured = as.character(slopes$meassured)

ggplot(slopes, aes(x=age, y=slope, colour=group))+
    # geom_point()+
    geom_jitter(alpha=0.4)+
    geom_smooth()+
    theme_bw()+
    scale_color_manual(values=as.character(sort(unique(slopes$colour))))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

tmp = subset(slopes, meassured == 100 | meassured == 200 | meassured == 300 | meassured == 400)
(slopesPlot = ggplot(tmp, aes(x=group, y=slope, colour=meassured))+
		geom_jitter(alpha=0.2)+
		geom_boxplot(alpha=0.6)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))

pdf(paste0(files$out, "_slopes_timeIntervals.pdf"), width=11.69, height=8.27, paper='special'); plot(slopesPlot); dev.off()

#----
#---- Get root dates for the different suergroups --------------------------------------------------

tmp = subset(data, N==1)
tmp = tmp[tmp$tree %in% subset(data %>% group_by(tree, root) %>% summarise(max=max(N)), max>300)$tree,]
tmp = tmp %>% group_by(tree, root) %>% summarise(min=min(time),
                                                 median=quantile(time, 0.50),
                                                 max=max(time))
(tmp = tmp[order(tmp$median, decreasing=TRUE),])

#----
#---- Summarise, subset and export table -----------------------------------------------------------

# Summarise table
data$time = data$time * files$factor
summ = data %>% group_by(tree, root, N, logN) %>% summarise(q05=quantile(time, 0.05),
                                                            median=quantile(time, 0.50),
                                                            q95=quantile(time, 0.95),
                                                            replicates=length(unique(file)))

# only use lineages and time points with at least 50% of the replicates
tmp = max(summ[!grepl("MC01", summ$tree),]$replicates)
summ = subset(summ, replicates >= tmp*0.5)

# Only use groups with at least 300 lineages
summ = summ[summ$tree %in% names(table(summ$tree)[table(summ$tree) > 300]),]

# Exporting the dataset ____________________________________________________________________________
write.table(summ, paste0(files$out, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Now get origin times of all trees ________________________________________________________________
tmp = subset(data, N==1)
write.table(tmp, paste0(files$out, "_LCA.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

tmp = subset(tmp, tree=="Main_tree")
round(tmp %>% summarise(median=median(time), min=min(time), max=max(time)))

#----
#---- Reading table if already available -----------------------------------------------------------

summ = fread(paste0(files$out, ".tsv"))

#----
#---- beautifying the dataset for plotting ---------------------------------------------------------

sort(unique(summ$tree))
summ$tree = factor(summ$tree, 
                   levels=c("Main_tree", "MC01-main_tree",
                            "Amoebozoa", "Nucletmycea", "Holozoa", 
                            "Metamonada", "Discoba", 
                            "Haptista", "Cryptista", "Archaeplastida", 
                            "Rhizaria", "Stramenopila", "Alveolata"))

summ$root = factor(summ$root, levels=c("rootA", "rootD"))

{summ$colour = as.character(summ$tree)
	summ$colour[which(summ$colour=="Main_tree")]="black"
	summ$colour[which(summ$colour=="MC01-main_tree")]="grey80"
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

#----
#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table
geos = subset(geo, time > min(-summ$q95, na.rm=TRUE))

(lttplot = ggplot(summ)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_line(aes(x=-median, y=N, color=tree)) +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    # geom_line(aes(x=-q05, y=logN, color=tree), alpha=0.2) +
    # geom_line(aes(x=-q95, y=logN, color=tree), alpha=0.2) +
    # geom_text(data=geos, aes(x=mid, y=max(summ$N)+10, label=era), size=2.5)+
    facet_wrap(~root, ncol=2)+
    geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    # scale_y_continuous(breaks=1:(round(max(summ$logN), 0))) +
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))
pdf(paste0(files$out, ".pdf"), width=11.69*1.5, height=8.27, paper='special'); plot(lttplot); dev.off()


(lttplotA = ggplot(subset(summ, root=="rootA"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_line(aes(x=-median, y=N, color=tree)) +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    # geom_line(aes(x=-q05, y=logN, color=tree), alpha=0.2) +
    # geom_line(aes(x=-q95, y=logN, color=tree), alpha=0.2) +
    # geom_text(data=geos, aes(x=mid, y=max(summ$N)+10, label=era), size=2.5)+
    facet_wrap(~root, ncol=2)+
    geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    # scale_y_continuous(breaks=1:(round(max(summ$logN), 0))) +
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))
pdf(paste0(files$out, "_rootA.pdf"), width=11.69, height=8.27, paper='special'); plot(lttplotA); dev.off()


(lttplotD = ggplot(subset(summ, root=="rootD"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_line(aes(x=-median, y=N, color=tree)) +
    geom_ribbon(aes(xmin=-q95, xmax=-q05, y=N, fill=tree), alpha=0.2, colour = NA) +
    # geom_line(aes(x=-q05, y=logN, color=tree), alpha=0.2) +
    # geom_line(aes(x=-q95, y=logN, color=tree), alpha=0.2) +
    # geom_text(data=geos, aes(x=mid, y=max(summ$N)+10, label=era), size=2.5)+
    facet_wrap(~root, ncol=2)+
    geom_text(data=geos, aes(x=mid, y=max(summ$N)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                       minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    # scale_y_continuous(breaks=1:(round(max(summ$logN), 0))) +
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)"))
pdf(paste0(files$out, "_rootD.pdf"), width=11.69, height=8.27, paper='special'); plot(lttplotD); dev.off()



#----
