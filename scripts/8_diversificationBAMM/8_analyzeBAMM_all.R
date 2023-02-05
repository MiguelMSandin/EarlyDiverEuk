#---- 
#---- loading packages ----

library(data.table)
library(dplyr)
library(ape)
library(phangorn)
library(ggplot2)

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid <- apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#---- 
# 'setwd' to the root directory where you have RTT and shift files
setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/")
#---- 
#---- Set file names -------------------------------------------------------------------------------

files = list(files=grep("rootD.*\\/bamm\\/.*RTT\\.tsv$", dir(recursive=TRUE), value=TRUE),
			 shifts=grep("rootD.*\\/bamm\\/.*shifts\\.tsv$", dir(recursive=TRUE), value=TRUE), 
			 trees=grep("rootD.*\\/clades\\/clade_.*\\.tre$", dir(recursive=TRUE), value=TRUE),
			 dirPlots="plots/",
			 RTTplot="bamm_RTT_allMedian.pdf",
			 shiftsPlot="bamm_shifts_all_violin.pdf")

if(!dir.exists(files$dirPlots)){dir.create(files$dirPlots)}

#---- 
#---- plot RTT BAMM output -------------------------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])

# Reading all RTT files
data = data.frame()
for(file in files$files){
	cat("\r  Working on file ", file, " (", grep(file, files$files), "/", length(files$files), ")          ", sep="")
	tmp = fread(file)
	tmp$file = file
	tmp$clade = file %>% sub(".*diversification_", "", .) %>% sub("_.*", "", .)
	tmp$diversity = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
	data = rbind(data, tmp)
};rm(file, tmp); cat("\n")

# Create another variable for minimum and maximum diversity
tmp = c()
for(cl in unique(data$clade)){
	ss = subset(data, clade==cl)
	tmp = c(tmp, ifelse(ss$diversity==min(ss$diversity), "min", "max"))
}; rm(cl, ss)
data$div = tmp; rm(tmp)

# Create a variable for the time point sampled per file
data$timePoint = rep(seq(1:100), length(unique(data$file))*length(unique(data$rate)))

# Summarise median slope per clade
summ = data %>% group_by(clade, div, rate, timePoint) %>% summarise(timeMedian=median(time),
																	timeSD=sd(time),
																	median2=median(median),
																	HPDlowMedian=median(HPDlow),
																	HPDuppMedian=median(HPDupp))

# Adding colours
{summ$colour = summ$clade
	summ$colour[which(summ$colour=="Amoebozoa")] =      "royalblue1"
	summ$colour[which(summ$colour=="Nucletmycea")] =    "steelblue2"
	summ$colour[which(summ$colour=="Holozoa")] =        "steelblue4"
	summ$colour[which(summ$colour=="Metamonada")] =     "forestgreen"
	summ$colour[which(summ$colour=="Discoba")] =        "orange2"
	summ$colour[which(summ$colour=="Haptista")] =       "yellow1"
	summ$colour[which(summ$colour=="Cryptista")] =      "hotpink1"
	summ$colour[which(summ$colour=="Archaeplastida")] = "darkseagreen3"
	summ$colour[which(summ$colour=="Rhizaria")] =       "darkorchid2"
	summ$colour[which(summ$colour=="Stramenopila")] =   "darkorchid3"
	summ$colour[which(summ$colour=="Alveolata")] =      "darkorchid4"}

# Adding factors for plotting
summ$rate = factor(summ$rate, levels=c("extinction","speciation","diversification"))
summ$clade = factor(summ$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
summ$colour = factor(summ$colour, levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2","darkorchid3","darkorchid4"))

# Plotting
geos <- subset(geo, time > -max(summ$timeMedian))

(rttPlot = ggplot(summ, aes(x=-timeMedian, y=median2, colour=clade))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_segment(aes(x=-timeMedian-timeSD/2, xend=-timeMedian+timeSD/2, y=median2, yend=median2), 
					 color="grey90", size=2, lineend="round")+
        geom_ribbon(aes(ymin=HPDlowMedian, ymax=HPDuppMedian, fill=colour, color=NULL), alpha=0.1)+
        geom_line()+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
		facet_grid(rate ~ div, scales="free")+
		scale_color_manual(values=as.character(sort(unique(summ$colour)))) +
		scale_fill_manual(values=as.character(sort(unique(summ$colour)))) +
        scale_x_continuous(breaks=seq(-max(round(summ$timeMedian, -2)), 0, 100),
                           minor_breaks=seq(-max(round(summ$timeMedian, -2)), 0, 50)) +
        labs(y="Rates Through Time", x="Time before present")+
		theme_classic())

# Export plot
pdf(paste0(files$dirPlots, "/", files$RTTplot), width=11.69, height=8.27, paper='special'); plot(rttPlot); dev.off()

#

#----  
#---- Analyze when diversification shifts tend to happen -------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files")])

files$trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE)
files

# Read datasets
data = data.frame()
for(file in files$shifts){
	cat("\r  Working on file ", file, " (", grep(file, files$shifts), "/", length(files$shifts), ")          ", sep="")
	tmp = fread(file)
	tmp$file = file
	tmp$clade = file %>% sub(".*diversification_", "", .) %>% sub("_div.*", "", .)
	tmp$diversity = as.numeric(file %>% sub(".*div", "", .) %>% sub("_.*", "", .))
	data = rbind(data, tmp)
}; rm(file, tmp); cat("\rDone", rep(" ", 100), "\n")

# Create another variable with minimum and maximum diversity




tmp = c()
for(clad in unique(data$clade)){
	ss = subset(data, clade==clad)
	m = mean(unique(ss$diversity))
	tmp = c(tmp, ifelse(ss$diversity < m, "min", "max"))
}; rm(clad, ss, m)
data$div = tmp; rm(tmp)





# Estimate the shifts on diversification rate
# lam(t) = lam1 * exp(lam2 * t)
# mu(t) = mu1 * exp(mu2 * t)
data$diversification = data$lam1 - data$mu2

# Get whether the shift corresponds to a decay or a growth
change = c()
for(f in unique(data$file)){
	cat("\r  Working on (", grep(f, unique(data$file)), "/", length(unique(data$file)), ") ", f, "                    ", sep="")
	dirRoot = f %>% sub("\\/[^\\/]+$", "", .) %>% sub("bamm/.*", "", .)
	clade = f %>% sub(".*diversification_", "", .) %>% sub("_div.*", "", .)
	ss = subset(data, file==f)
	tree = files$trees %>% grep(dirRoot, ., value=TRUE) %>% grep(clade, ., value=TRUE)
	tree = read.tree(tree)
	changeF = c("root")
	for(n in 2:length(ss$node)){
		noden = ss$node[n]
		dn = subset(ss, node==noden)$diversification
		ancestors = Ancestors(tree, n)
		parent = max(ancestors[ancestors %in% ss$node])
		dp = subset(ss, node==parent)$diversification
		tmp = ifelse(dn > dp, "growth", "decay")
		changeF = c(changeF, tmp)
	}
	change = c(change, changeF)
}; rm(f, dirRoot, clade, ss, tree, changeF, n, noden, dn, ancestors, parent, dp, tmp); cat("\rDone", rep(" ", 100), "\n")

data$change = change; rm(change)
data$change = factor(data$change, levels=c("root", "growth", "decay"))

# Remove shifts happening in 4 or less phylogenetic trees during a 10 years span
# This is done to compensate the low ESS for big clades (e.g.: Holozoa, Nucletmycea and Alveolata)
span = 50
data$replicates = 0
for(c in unique(data$clade)){
	for(d in unique(data$div)){
		ss = subset(data, clade == c & div == d)
		ss = subset(ss, change!="root")
		sequence = seq(min(ss$timeRev), 0, span)
		for(s in 1:(length(sequence)-1)){
			cat("\r  Working on ", c, ": Analyzing '", d, "' diversity: ", round(s/length(sequence)*100), "%", sep="")
			sss = subset(ss, timeRev >= sequence[s] & timeRev <= sequence[s+1])
			if(nrow(sss) > 0){
				for(n in 1:nrow(sss)){
					data$replicates[which(data$node==sss$node[n] & data$file==sss$file[n])] = length(unique(sss$file))
				}
			}
		}
	}; cat("\r  Working on ", c, ": Done                              \n", sep="")
}; rm(c, d, ss, sequence, s, sss, n)


table(data$replicates)




# And plot
geos <- subset(geo, time>-max(data$time))

# Add factors
data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

(shifts = ggplot(data, aes(x=timeRev, y=clade, colour=change))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_violin(scale="width")+
		facet_wrap(~div, ncol=1)+
		scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
						   minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
		labs(title="Shifts in diversification rates",
			 y="Clades", x="Time before present (Ma)")+
		scale_colour_manual(values=c("grey20", "springgreen4", "orangered4")) +
		theme_classic())

# And export
pdf(paste0(files$dirPlots, "/", files$shiftsPlot), width=11.69, height=8.27, paper='special'); plot(shifts); dev.off()

# And now only shifts happening in the proterozoic
geos <- subset(geo, time < -540 & time>-max(data$time))

(shiftsP = ggplot(subset(data, timeRev < -541), aes(x=timeRev, y=clade, fill=change))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_violin(scale="width", alpha=0.4)+
		facet_wrap(~div, ncol=1)+
		scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
						   minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
		labs(title="Shifts in diversification rates on the Proterozoic",
			 y="Clades", x="Time before present (Ma)")+
		scale_fill_manual(values=c("grey20", "springgreen4", "orangered4")) +
		theme_classic())

pdf(paste0(files$dirPlots, "/", sub("\\.pdf", "_proterozoic.pdf", files$shiftsPlot)), width=11.69, height=8.27, paper='special'); plot(shiftsP); dev.off()

#

#----  