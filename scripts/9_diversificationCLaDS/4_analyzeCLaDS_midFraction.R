#---- 
#---- loading packages ----

library(HDInterval)
library(data.table)
library(dplyr)
library(ggplot2)
 
geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----  
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")
rm(list=ls()[!ls() %in% c("geo")])
# .rs.restartR()
#---- Set file names -------------------------------------------------------------------------------

# files = list(folder="rootA/fRAng2-2-1",
#              outPrefix="plots/cladsMid/divMid_fRAng2-2-1")
files = list(folder="rootD/fRAgamma2-1",
             outPrefix="plots/cladsMid/divMid_fRAgamma2-1")

{files$files=grep(paste0(files$folder, ".*clads.*RTT\\.tsv"), dir(recursive=TRUE), value=TRUE)
  files$param=grep(paste0(files$folder, ".*clads.*parameters.txt"), dir(recursive=TRUE), value=TRUE)
  files$tipRates=grep(paste0(files$folder, ".*clads.*tipRates.txt"), dir(recursive=TRUE), value=TRUE)
  files$outDir = dirname(files$outPrefix)
  files}

if(!dir.exists(files$outDir)){dir.create(files$outDir)}

#----
#---- Read all model parameters files --------------------------------------------------------------

# Sigma (σ): Stochastic parameters
# Alpha (α): trend at speciation (whether daughter rates tend to be higher or lower than parental rates)
# Epsilon (ε): Extinction rates constant turn-over (that is, μi/λi = ε for all lineages; ClaDS2)
# Lambda0 (λ0): Speciation rate at the beginning of the process
# mu0 (μ0): Extinction rate at the beginning of the process

# Now extract all files
param = data.table()
for(file in files$param){
  cat("\r  Loading file (", grep(file, files$param), "/", length(files$param), ") ", file, "                    " , sep="", end="")
  tmp = fread(file)
  tmp$file = file
  tmp$logMeanLambdaTip = log10(tmp$meanLambdaTip)
  tmp$m = tmp$alpha * exp((tmp$sigma^2)/2)
  tmp$root = sub("\\/(f|r)RA.*", "", file)
  tmp$replicate = file %>% sub("root.\\/", "", .) %>% sub("\\/clads.*", "", .)
  tmp$group = file %>% sub(".*clads(Mid)*\\/", "", .) %>% sub("\\/.*", "", .)
  tmp$div = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
  param = rbind(param, tmp)
}; rm(tmp); cat("\n")

tiprates = data.table()
for(file in files$tipRates){
  cat("\r  Loading file (", grep(file, files$tipRates), "/", length(files$tipRates), ") ", file, "                    " , sep="", end="")
  tmp = fread(file)
  tiprates = rbind(tiprates, data.frame(root=sub("\\/(f|r)RA.*", "", file),
                                        replicate = file %>% sub("root.\\/", "", .) %>% sub("\\/clads.*", "", .),
                                        group = file %>% sub(".*clads(Mid)*\\/", "", .) %>% sub("\\/.*", "", .),
                                        div = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .)),
                                        medianLambdaTip=median(tmp$V1)))
}; rm(tmp); cat("\n")

param = merge(param, tiprates, by=c("root", "replicate", "group", "div"))

param$diversity = NA
for(cl in unique(param$group)){
	div = subset(param, group==cl)$div
	m = round(mean(div))
	for(d in div){
	  tmp = ifelse(d > m, "max", ifelse(d < m, "min", "mean"))
	  param$diversity[which(param$group == cl & param$div == d)] = tmp
	}
};rm(cl, div, m, d, tmp)

# param$mu0 = param$epsilon * param$lambda0
# param$m = param$alpha * exp((param$sigma^2) / 2)

#---- Prepare for plotting -------------------------------------------------------------------------
paramm = melt(param, id.vars=c("file", "root", "replicate", "group", "div", "diversity"))

{paramm$colour = paramm$group
	paramm$colour[which(paramm$colour=="Amoebozoa")]="royalblue1"
	paramm$colour[which(paramm$colour=="Nucletmycea")]="steelblue2"
	paramm$colour[which(paramm$colour=="Holozoa")]="steelblue4"
	paramm$colour[which(paramm$colour=="Metamonada")]="forestgreen"
	paramm$colour[which(paramm$colour=="Discoba")]="orange2"
	paramm$colour[which(paramm$colour=="Haptista")]="yellow1"
	paramm$colour[which(paramm$colour=="Cryptista")]="hotpink1"
	paramm$colour[which(paramm$colour=="Archaeplastida")]="darkseagreen3"
	paramm$colour[which(paramm$colour=="Rhizaria")]="darkorchid2"
	paramm$colour[which(paramm$colour=="Stramenopila")]="darkorchid3"
	paramm$colour[which(paramm$colour=="Alveolata")]="darkorchid4"}

paramm$group = factor(paramm$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

paramm$colour = factor(paramm$colour,
					levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))

paramm$diversity = factor(paramm$diversity, levels=c("max", "mean","min"))

paramm$grouping = paste0(paramm$group, "-", paramm$diversity)

#
#---- Plotting -------------------------------------------------------------------------------------

paramms = subset(paramm, variable!="meanLambdaTip" & variable!="logMeanLambdaTip")

(paramDens = ggplot(paramms, aes(x=value, y=group, fill=diversity, colour=diversity))+
   geom_point(alpha=0.6)+
   facet_wrap( ~ variable, scales="free_x", axes="all_x",
              nrow=length(unique(paramms$group)), ncol=length(unique(paramms$variable)))+
   scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107")) +
   scale_colour_manual(values=c("#D81B60", "#1E88E5", "#FFC107")) +
   theme_minimal())

# (paramDens = ggplot(paramm, aes(x=value, fill=group, colour=group))+
#  	geom_density(alpha=0.6)+
# 	# geom_freqpoly(alpha=0.6, position="identity")+
# 	facet_wrap(variable ~ diversity, scales="free", ncol=length(unique(paramm$diversity)), nrow=length(unique(paramm$variable)))+
# 	scale_fill_manual(values=as.character(sort(unique(paramm$colour)))) +
# 	scale_colour_manual(values=as.character(sort(unique(paramm$colour)))) +
# 	theme_minimal())

pdf(paste0(files$outPrefix, "_parameters.pdf"), width=11.69, height=11.69, paper='special'); plot(paramDens); dev.off()

#----
#---- Read all RTT and DTT extracted files from the Rdata files ------------------------------------

# Get whether the files correspond to the maximum, mean or minimum diversity
files$clades = unique(files$files %>% sub(".*clads(Mid)*\\/", "", .) %>% sub("\\/.*", "", .))
files$diversity = data.frame()
for(clade in files$clades){
	f = grep(clade, files$files, value=TRUE)
	div = f %>% sub(".*_div", "", .) %>% sub("_.*", "", .)
	div = as.numeric(div)
	d = data.frame(file=f, div=ifelse(div == max(div), "max", ifelse(div == min(div), "min", "mean")))
	files$diversity = rbind(files$diversity, d)
};rm(clade, f, div, d)

files

# Now extract all files
rtt = data.frame()
dtt = data.frame()
for(file in files$files){
    cat("\r  Loading file (", grep(file, files$files), "/", length(files$files), ") ", file, "                    " , sep="", end="")
    tmp = fread(file)
    div = files$diversity[grep(file, files$diversity[,1]),2]
    tmp$file = file
    tmp$root = sub("\\/(f|r)RA.*", "", file)
    tmp$group = file %>% sub(".*clads(Mid)*\\/", "", .) %>% sub("\\/.*", "", .)
    tmp$revTime = -rev(tmp$time)
    tmp$div = div
    tmp$timePoint = 1:nrow(tmp)
    rtt = rbind(rtt, tmp)
    
    tmp = fread(gsub("RTT", "DTT", file))
    tmp$file = file
    tmp$root = sub("\\/(f|r)RA.*", "", file)
    tmp$group = file %>% sub(".*clads(Mid)*\\/", "", .) %>% sub("\\/.*", "", .)
    tmp$div = div
    tmp$timePoint = 1:nrow(tmp)
    dtt = rbind(dtt, tmp)
}; rm(tmp, div); cat("\n")

#---- Prepare for plotting -------------------------------------------------------------------------

geos = subset(geo, time>min(rtt$revTime))

{rtt$colour = rtt$group
    rtt$colour[which(rtt$colour=="Amoebozoa")]="royalblue1"
    rtt$colour[which(rtt$colour=="Nucletmycea")]="steelblue2"
    rtt$colour[which(rtt$colour=="Holozoa")]="steelblue4"
    rtt$colour[which(rtt$colour=="Metamonada")]="forestgreen"
    rtt$colour[which(rtt$colour=="Discoba")]="orange2"
    rtt$colour[which(rtt$colour=="Haptista")]="yellow1"
    rtt$colour[which(rtt$colour=="Cryptista")]="hotpink1"
    rtt$colour[which(rtt$colour=="Archaeplastida")]="darkseagreen3"
    rtt$colour[which(rtt$colour=="Rhizaria")]="darkorchid2"
    rtt$colour[which(rtt$colour=="Stramenopila")]="darkorchid3"
    rtt$colour[which(rtt$colour=="Alveolata")]="darkorchid4"}

rtt$group = factor(rtt$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

rtt$colour = factor(rtt$colour,
                     levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))

# rtt$root = factor(rtt$root, levels=c("rootD", "rootA"))
rtt$div = factor(rtt$div, levels=c("max", "mean", "min"))
rtt$grouping = paste0(rtt$group, "-", rtt$div)

{dtt$colour = dtt$group
	dtt$colour[which(dtt$colour=="Amoebozoa")]="royalblue1"
	dtt$colour[which(dtt$colour=="Nucletmycea")]="steelblue2"
	dtt$colour[which(dtt$colour=="Holozoa")]="steelblue4"
	dtt$colour[which(dtt$colour=="Metamonada")]="forestgreen"
	dtt$colour[which(dtt$colour=="Discoba")]="orange2"
	dtt$colour[which(dtt$colour=="Haptista")]="yellow1"
	dtt$colour[which(dtt$colour=="Cryptista")]="hotpink1"
	dtt$colour[which(dtt$colour=="Archaeplastida")]="darkseagreen3"
	dtt$colour[which(dtt$colour=="Rhizaria")]="darkorchid2"
	dtt$colour[which(dtt$colour=="Stramenopila")]="darkorchid3"
	dtt$colour[which(dtt$colour=="Alveolata")]="darkorchid4"}

dtt$group = factor(dtt$group, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

dtt$colour = factor(dtt$colour,
					 levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2","darkorchid3","darkorchid4"))

# dtt$root = factor(dtt$root, levels=c("rootD", "rootA"))
dtt$div = factor(dtt$div, levels=c("max", "mean", "min"))
dtt$grouping = paste0(dtt$group, "-", dtt$div)

# write.table(rtt, paste0(files$outDir, "RTT", sub("pdf$", "tsv", files$outSuffix)), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
# write.table(dtt, paste0(files$outDir, "DTT", sub("pdf$", "tsv", files$outSuffix)), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#
#---- Plotting RTT ---------------------------------------------------------------------------------

(rttplot = ggplot(rtt)+
   geom_vline(xintercept = geos$time, color="lightgrey") +
   geom_ribbon(aes(x=revTime, ymin=HPD_05, ymax=HPD_95, fill=colour, color=NULL, group=grouping), alpha=0.05)+
   # geom_line(aes(x=revTime, y=rate, colour=group, group=grouping, linetype=div))+
   geom_line(aes(x=revTime, y=rate, colour=group, group=grouping, linewidth=div))+
   scale_y_log10() + annotation_logticks(sides = 'l')+
   scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
   scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
   # scale_linetype_manual(values=c("dashed", "solid", "dotted")) +
   scale_linewidth_manual(values=c(0.2, 1, 0.4)) +
   scale_x_continuous(breaks=seq(min(round(rtt$revTime, -2)), 0, 100),
                      minor_breaks=seq(min(round(rtt$revTime, -2)), 0, 50)) +
   labs(y="Rate Through Time (Ln)", x="Time (Ma)")+
   theme_classic()+
   theme(legend.position = "none"))

pdf(paste0(files$outPrefix, "_RTT.pdf"), width=11.69, height=8.27, paper='special'); plot(rttplot); dev.off()

#
#---- Plotting DTT ---------------------------------------------------------------------------------

(dttplot = ggplot(dtt)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_line(aes(x=time, y=lineages, colour=group, group=grouping, linewidth=div))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    scale_color_manual(values=as.character(sort(unique(dtt$colour)))) +
    scale_fill_manual(values=as.character(sort(unique(dtt$colour)))) +
    scale_x_continuous(breaks=seq(min(round(dtt$time, -2)), 0, 100),
                       minor_breaks=seq(min(round(dtt$time, -2)), 0, 50)) +
    scale_linewidth_manual(values=c(0.2, 1, 0.4)) +
    labs(y="Estimated diversity through time (Ln)", x="Time (Ma)")+
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$outPrefix, "_DTT.pdf"), width=11.69, height=8.27, paper='special'); plot(dttplot); dev.off()

#
