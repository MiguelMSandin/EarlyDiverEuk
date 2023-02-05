#! /bin/Rscript

# loading packages _________________________________________________________________________________

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(RPANDA))
suppressPackageStartupMessages(library(HDInterval))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, 
						  -1600 , -1400, -1200, 
						  -1000, -720, -635, 
						  -541, -485.4, -443.8, -419, -358.9, -298.9, 
						  -251.9, -201.4, -145, 
						  -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", 
				  		 "Calymnian","Ectasian","Stenian",
				  		 "Tonian","Cryogenian","Ediacran",
				  		 "Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian",
				  		 "Triassic","Jurassic","Cretaceous",
				  		 "Paleogene","Neogene","Quaternary"),
				  era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic",
				  	  "Mesoproterozoic","Mesoproterozoic","Mesoproterozoic",
				  	  "Neoproterozoic","Neoproterozoic","Neoproterozoic",
				  	  "Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
				  	  "Mesozoic","Mesozoic","Mesozoic",
				  	  "Cenozoic","Cenozoic","Cenozoic"))

# Set parser _______________________________________________________________________________________
parser <- OptionParser()

parser <- add_option(parser, c("-f", "--rdata"), dest="file", type="character",
                     help="RData file name saved from CLaDS function in julia.")

parser <- add_option(parser, c("-b", "--burnin"), dest="burnin", type="integer", default=25,
					 help="The burnin (in percentage; 0-100), default='25'.")

parser <- add_option(parser, c("-o", "--out"), dest="out", type="character", default="NONE",
					 help="The output data prefix. By default, will remove the extension to the input file.")

parser <- add_option(parser, c("-v", "--verbose"), dest="verbose", action="store_true", default=TRUE,
                     help="If selected, will not print information to the console.")

args = parse_args(parser)

# Setting file names _______________________________________________________________________________
if(args$out== "NONE"){
	prefix = sub("\\.[^\\.]+$", "", args$file)
}else{
	prefix = args$out
}
outTree = paste0(prefix, "_RTTtree.pdf")
outRTTplot = paste0(prefix, "_RTT.pdf")
outRTTtable = paste0(prefix, "_RTT.tsv")
outLTTplot = paste0(prefix, "_DTT.pdf")
outLTTtable = paste0(prefix, "_DTT.tsv")

# Open the analyzed tree ___________________________________________________________________________
if(args$verbose){cat("  Reading file:", args$file, "\n")}
load(args$file)


# Checking output from ClaDS Julia _________________________________________________________________
if(args$verbose){cat("  Plotting diversification rates in the tree\n")}

tree = CladsOutput$tree
rates = CladsOutput$lambdai_map

pdf(outTree, width=11.69, height=8.27, paper='special')
plot_ClaDS_phylo(tree, rates, show.tip.label=TRUE, cex=0.2)
title(prefix)
tmp = dev.off()

# Plotting the rate through time plot ______________________________________________________________
if(args$verbose){cat("  Plotting diversification Rate Through Time plot\n")}
time <- apply(data.frame(CladsOutput$time_points[-length(CladsOutput$time_points)], CladsOutput$time_points[-1]),1 , mean)
clads_rtt <- data.frame(time=time, 
						rate=CladsOutput$RTT_map)

chains <- data.frame(time=time)
burning = length(CladsOutput$rtt_chains[[1]])*(args$burnin/100)
c <- 0
for(chain in CladsOutput$rtt_chains){
	c <- c + 1
	chaini <- as.data.frame(chain)
	colnames(chaini) <- paste0("Chain", c, "_iter", 1:ncol(chaini))
	chaini <- chaini[,-c(1:burning)]
	chains <- cbind(chains, chaini)
}; rm(c, chaini)
chains$toRemove <- NULL

hpd <- as.data.frame(t(apply(t(chains), 2, hdi)))
clads_rtt$HPD_05 <- hpd$lower
clads_rtt$HPD_95 <- hpd$upper
clads_rtt$mean = apply(t(chains), 2, mean)

geos <- subset(geo, time>-max(clads_rtt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

rttPlot <-ggplot(clads_rtt, aes(x=-rev(time), y=rate))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_line(linewidth=2)+
		geom_line(aes(y=HPD_05), colour="lightblue")+
		geom_line(aes(y=HPD_95), colour="lightblue")+
		# geom_line(aes(y=mean), colour="springgreen2")+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq((round(max(clads_rtt$time)*-100, -2)), 0, 100),
						   minor_breaks=seq((round(min(clads_rtt$time)*-100, -2)), 0, 50)) +
		theme_classic()+
		labs(title=paste0(prefix),
			 y="Diversification rate (Ln)", x="Time (Ma)")

pdf(outRTTplot, width=11.69, height=8.27, paper='special')
plot(rttPlot)
tmp = dev.off()
write.table(clads_rtt, outRTTtable, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Plotting the Lineages through time plot __________________________________________________________
if(args$verbose){cat("  Plotting Lineages Through Time plot\n")}

clads_ltt <- data.frame(time_points=CladsOutput$time_points,
						time=-rev(CladsOutput$time_points),
						lineages=CladsOutput$DTT_mean)

enhanced = data.frame(time=-rev(CladsOutput$time_points))
for(t in 1:length(CladsOutput$enhanced_trees)){
	tmp = as.data.frame(ltt.plot.coords(CladsOutput$enhanced_trees[[t]]$tree))
	tmp1 = c()
	for(ti in enhanced$time){
		i = ifelse(ti < min(tmp$time), 1, max(subset(tmp, time < ti)$N)+1)
		tmp1 = c(tmp1, i)
	}
	enhanced = cbind(enhanced, tmp1)
}; rm(t, tmp, tmp1, ti); cat("\n")
colnames(enhanced) = c("time", paste0("tree", 1:(ncol(enhanced)-1)))

clads_ltt$mean = apply(enhanced[,-1], 1, mean)
clads_ltt$min = apply(enhanced[,-1], 1, min)
clads_ltt$max = apply(enhanced[,-1], 1, max)

geos <- subset(geo, time>min(clads_ltt$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

lttPlot <-ggplot(clads_ltt, aes(x=time, y=lineages))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		# geom_point()+
		# geom_line(aes(y=min), colour="lightblue")+
		# geom_line(aes(y=max), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree1), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree2), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree3), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree4), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree5), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree6), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree7), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree8), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree9), colour="lightblue")+
		geom_line(data=enhanced, aes(x=time, y=tree10), colour="lightblue")+
		geom_line(linewidth=2)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq((round(min(clads_ltt$time), -2)), 0, 100),
						   minor_breaks=seq((round(min(clads_ltt$time), -2)), 0, 50)) +
		theme_classic()+
		labs(title=paste0(prefix),
			 y="Diversity Through Time (Ln)", x="Time (Ma)")

pdf(outLTTplot, width=11.69, height=8.27, paper='special')
plot(lttPlot)
tmp = dev.off()
write.table(clads_ltt, outLTTtable, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Finished _________________________________________________________________________________________
if(args$verbose){cat("Done\n")}
