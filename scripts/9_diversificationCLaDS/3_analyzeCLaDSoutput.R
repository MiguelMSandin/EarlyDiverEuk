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
setwd("~/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/")
rm(list=ls()[!ls() %in% c("geo")])
# .rs.restartR()
#---- Set file names -------------------------------------------------------------------------------

files = list(files=grep("root.\\/.*clads\\/.*RTT\\.tsv", dir(recursive=TRUE), value=TRUE),
             param=grep("root.\\/.*clads\\/.*parameters.txt", dir(recursive=TRUE), value=TRUE),
             tipRates=grep("root.\\/.*clads\\/.*tipRates.txt", dir(recursive=TRUE), value=TRUE),
             outDir="plots/clads/",
             outPreffix="MC01-MC02_ClaDS_")

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
  tmp$calibration = sub("\\/root.*", "", file)
  tmp$root = file %>% sub("\\/(f|r)RA.*", "", .) %>% sub(".*\\/", "", .)
  tmp$replicate = file %>% sub(".*root.\\/", "", .) %>% sub("\\/clads.*", "", .)
  tmp$group = file %>% sub(".*clads\\/", "", .) %>% sub("\\/.*", "", .)
  tmp$div = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
  param = rbind(param, tmp)
}; rm(tmp); cat("\n")

tiprates = data.table()
for(file in files$tipRates){
  cat("\r  Loading file (", grep(file, files$tipRates), "/", length(files$tipRates), ") ", file, "                    " , sep="", end="")
  tmp = fread(file)
  tiprates = rbind(tiprates, data.frame(calibration = sub("\\/root.*", "", file),
                                        root = file %>% sub("\\/(f|r)RA.*", "", .) %>% sub(".*\\/", "", .),
                                        replicate = file %>% sub(".*root.\\/", "", .) %>% sub("\\/clads.*", "", .),
                                        group = file %>% sub(".*clads\\/", "", .) %>% sub("\\/.*", "", .),
                                        div = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .)),
                                        medianLambdaTip=median(tmp$V1)))
}; rm(tmp); cat("\n")

param = merge(param, tiprates, by=c("calibration", "root", "replicate", "group", "div"))

param$diversity = NA
for(cl in unique(param$group)){
	ss = subset(param, group==cl)
	tmp = mean(ss$div)
	for(d in unique(ss$div)){
		param$diversity[which(param$group == cl & param$div == d)] = ifelse(d < tmp, "min", "max")
	}
};rm(cl, ss, tmp, d)

# param$mu0 = param$epsilon * param$lambda0
# param$m = param$alpha * exp((param$sigma^2) / 2)

table(param$replicate)
param$calibration = sub("-2", "", param$calibration)


#---- Prepare for plotting -------------------------------------------------------------------------
paramm = melt(param, id.vars=c("file", "calibration", "root", "replicate", "group", "div", "diversity"))

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

paramm$calibration = factor(paramm$calibration, levels=c("MC01", "MC02"))
paramm$root = factor(paramm$root, levels=c("rootD", "rootA"))

paramm$wrapping = paste0(paramm$root, "-", paramm$diversity)
paramm$wrapping = factor(paramm$wrapping, levels=c("rootD-max", "rootD-min", "rootA-max", "rootA-min"))

paramm$grouping = paste0(paramm$group, "-", paramm$diversity)
paramm$grouping = factor(paramm$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
                                                   "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
                                                   "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
                                                   "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))

#
#---- Get some stats -------------------------------------------------------------------------------

paramm %>% group_by(variable) %>% summarise(median(value))

tmp = paramm %>% group_by(group, calibration, root, diversity, variable) %>% summarise(min=min(value), q05=quantile(value, 0.05), q25=quantile(value, 0.25),
                                                                          mean=mean(value), sd=sd(value),median=median(value),
                                                                          q75=quantile(value, 0.75), q95=quantile(value, 0.95), max=max(value))

tmp = paramm %>% group_by(group, calibration, variable) %>% summarise(min=min(value), q05=quantile(value, 0.05), q25=quantile(value, 0.25),
                                                         mean=mean(value), sd=sd(value),median=median(value),
                                                         q75=quantile(value, 0.75), q95=quantile(value, 0.95), max=max(value))

tmp1 = subset(tmp, variable=="medianLambdaTip")
tmp1 = subset(tmp, variable=="sigma")

#---- Plotting -------------------------------------------------------------------------------------

paramms = subset(paramm, variable!="meanLambdaTip" & variable!="logMeanLambdaTip")

(paramDens = ggplot(paramms, aes(x= value, y=group, fill=group, colour=group, group=grouping))+
    geom_violin(aes(alpha=grouping))+
    facet_wrap(root+calibration ~ variable, scales="free_x", axes="all_x",
                nrow=length(unique(c(paramms$diversity, param$root))), ncol=length(unique(paramms$variable)))+
    scale_fill_manual(values=as.character(sort(unique(paramms$colour)))) +
    scale_colour_manual(values=as.character(sort(unique(paramms$colour)))) +
    scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
    theme_minimal()+
    theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "parameters.pdf"), width=11.69, height=11.69, paper='special'); plot(paramDens); dev.off()

#----
#---- Read all RTT and DTT extracted files from the Rdata files ------------------------------------

# Get whether the files correspond to the maximum or minimum diversity
files$clades = unique(files$files %>% sub(".*clads\\/", "", .) %>% sub("\\/.*", "", .))
files$diversity = data.frame()
for(clade in files$clades){
	f = grep(clade, files$files, value=TRUE)
	div = f %>% sub(".*_div", "", .) %>% sub("_.*", "", .)
	div = as.numeric(div)
	m = mean(div)
	d = data.frame(f, fifelse(div > m, "max", "min"))
	files$diversity = rbind(files$diversity, d)
};rm(clade, f, div, m, d)

files

# Now extract all files
rtt_raw = data.frame()
dtt_raw = data.frame()
for(file in files$files){
    cat("\r  Loading file (", grep(file, files$files), "/", length(files$files), ") ", file, "                    " , sep="", end="")
    tmp = fread(file)
    div = files$diversity[grep(file, files$diversity[,1]),2]
    tmp$file = file
    tmp$calibration = file %>% sub("\\/root.*", "", .) %>% sub("-2", "", .)
    tmp$root = file %>% sub("\\/(f|r)RA.*", "", .) %>% sub(".*\\/", "", .)
    tmp$replicate = file %>% sub(".*root.\\/", "", .) %>% sub("\\/clads.*", "", .)
    tmp$group = file %>% sub(".*clads\\/", "", .) %>% sub("\\/.*", "", .)
    tmp$revTime = -rev(tmp$time)
    tmp$div = div
    tmp$timePoint = 1:nrow(tmp)
    rtt_raw = rbind(rtt_raw, tmp)
    
    tmp = fread(gsub("RTT", "DTT", file))
    tmp$file = file
    tmp$calibration = file %>% sub("\\/root.*", "", .) %>% sub("-2", "", .)
    tmp$root = file %>% sub("\\/(f|r)RA.*", "", .) %>% sub(".*\\/", "", .)
    tmp$replicate = file %>% sub(".*root.\\/", "", .) %>% sub("\\/clads.*", "", .)
    tmp$group = file %>% sub(".*clads\\/", "", .) %>% sub("\\/.*", "", .)
    tmp$div = div
    tmp$timePoint = 1:nrow(tmp)
    dtt_raw = rbind(dtt_raw, tmp)
}; rm(tmp, div); cat("\n")

#
#---- Summarise all slopes -------------------------------------------------------------------------

rtt = rtt_raw %>% group_by(calibration, group, root, div, timePoint) %>% 
  summarise(timeMedian=median(revTime),
            timeHPD05=hdi(revTime)[1],
            timeHPD95=hdi(revTime)[2],
            median2=median(rate),
            HPD05median=median(HPD_05),
            HPD95median=median(HPD_95))

dtt = dtt_raw %>% group_by(calibration, group, root, div, timePoint) %>% 
  summarise(timeMedian=median(time),
            timeHPD05=hdi(time)[1],
            timeHPD95=hdi(time)[2],
            median2=median(lineages),
            HPD05median=median(min),
            HPD95median=median(max))

#
#---- Prepare for plotting -------------------------------------------------------------------------

geos = subset(geo, time>min(rtt$timeHPD05))

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

rtt$calibration = factor(rtt$calibration, levels=c("MC01", "MC02"))
rtt$root = factor(rtt$root, levels=c("rootD", "rootA"))

rtt$grouping = paste0(rtt$group, "-", rtt$div)
unique(rtt$grouping)
rtt$grouping = factor(rtt$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
                                             "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
                                             "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
                                             "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))


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


dtt$calibration = factor(dtt$calibration, levels=c("MC01", "MC02"))
dtt$root = factor(dtt$root, levels=c("rootD", "rootA"))

dtt$grouping = paste0(dtt$group, "-", dtt$div)
unique(dtt$grouping)
dtt$grouping = factor(dtt$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
                                             "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
                                             "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
                                             "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))

write.table(rtt, paste0(files$outDir, "/", files$outPreffix, "RTT.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(dtt, paste0(files$outDir, "/", files$outPreffix, "DTT.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#
#---- Plotting RTT ---------------------------------------------------------------------------------

(rttplot = ggplot(rtt)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=timeMedian, yend=timeMedian),
                 color="grey90", linewidth=2, lineend="round")+
    geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, fill=colour, color=NULL, group=grouping), alpha=0.1)+
    geom_line(aes(x=timeMedian, y=median2, colour=group, group=grouping, alpha=grouping))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(calibration ~ root, scales="free")+
    scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
    scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
                       minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
    labs(y="Rate Through Time (Ln)", x="Time (Ma)")+
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "RTT.pdf"), width=8.27, height=8.27, paper='special'); plot(rttplot); dev.off()

#
#---- Plotting DTT ---------------------------------------------------------------------------------

(dttplot = ggplot(dtt)+
   geom_vline(xintercept = geos$time, color="lightgrey") +
   geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=median2, yend=median2), 
                color="grey90", linewidth=2, lineend="round")+
   geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, fill=colour, color=NULL, group=grouping), alpha=0.1)+
   geom_line(aes(x=timeMedian, y=median2, colour=group, group=grouping, alpha=grouping))+
   scale_y_log10() + annotation_logticks(sides = 'l')+
   facet_grid(calibration ~ root, scales="free")+
   scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
   scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
   scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
   scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
                      minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
   labs(y="Estimated diversity through time (Ln)", x="Time (Ma)")+
   theme_classic()+
   theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "DTT.pdf"), width=11.69, height=8.27, paper='special'); plot(dttplot); dev.off()

(dttplotm = ggplot(subset(dtt, div=="min"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=median2, yend=median2), 
                 color="grey90", linewidth=2, lineend="round")+
    geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, fill=colour, color=NULL, group=grouping), alpha=0.1)+
    geom_line(aes(x=timeMedian, y=median2, colour=group, group=grouping))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(calibration ~ root, scales="free")+
    scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
                       minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
    labs(y="Estimated diversity through time (Ln)", x="Time (Ma)")+
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "DTT_divMin.pdf"), width=11.69, height=8.27, paper='special'); plot(dttplotm); dev.off()

(dttplotM = ggplot(subset(dtt, div=="max"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=median2, yend=median2), 
                 color="grey90", linewidth=2, lineend="round")+
    geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, fill=colour, color=NULL, group=grouping), alpha=0.1)+
    geom_line(aes(x=timeMedian, y=median2, colour=group, group=grouping, alpha=0.6))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(calibration ~ root, scales="free")+
    scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
                       minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
    labs(y="Estimated diversity through time (Ln)", x="Time (Ma)")+
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "DTT_divMax.pdf"), width=11.69, height=8.27, paper='special'); plot(dttplotM); dev.off()

#
# #---- Plotting RTT and DTT -------------------------------------------------------------------------
# 
# tmp1 = rtt
# tmp1$data = "rtt"
# tmp2 = dtt
# tmp2$data = "dtt"
# tmp2$median2 = tmp2$median2-1
# tmp2$HPD05median = tmp2$HPD05median-1
# tmp2$HPD95median = tmp2$HPD95median-1
# 
# data = rbind(tmp1, tmp2); rm(tmp1, tmp2)
# data$data = factor(data$data, levels=c("rtt", "dtt"))
# 
# data$grouping = paste0(data$group, "-", data$div)
# unique(data$grouping)
# data$grouping = factor(data$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
#                                                "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
#                                                "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
#                                                "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))
# 
# (rttdttplot = ggplot(data)+
# 		geom_vline(xintercept = geos$time, color="lightgrey") +
# 		geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=median2, yend=median2), 
# 					 color="grey80", linewidth=2, lineend="round", alpha=0.1)+
# 		geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, fill=colour, color=NULL, group=grouping), alpha=0.1)+
# 		geom_line(aes(x=timeMedian, y=median2, colour=group, group=grouping, alpha=grouping))+
# 		scale_y_log10() + annotation_logticks(sides = 'l')+
# 		facet_grid(data ~ root, scales="free")+
# 		scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
# 		scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
#     scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
# 		scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
# 						   minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
# 		labs(x="Time (Ma)")+
# 		theme_classic()+
# 		theme(legend.position = "none"))
# 
# pdf(paste0(files$outDir, "/", files$outPreffix, "DTT-RTT.pdf"), width=11.69, height=11.69, paper='special'); plot(rttdttplot); dev.off()
# 
# # And plot fig 3 but for independent supergroups
# unique(data$group)
# 
# pdf(paste0(files$outDir, "RTT-DTT_independentSupergroups", files$outSuffix), width=11.69, height=8.27, paper='special')
# for(i in c("Discoba", "Metamonada", 
# 		   "Amoebozoa", "Nucletmycea", "Holozoa", 
# 		   "Archaeplastida", "Haptista", "Cryptista", "Rhizaria", "Alveolata", "Stramenopila")){
# 	tmp = subset(data, group==i)
# 	fig3i = ggplot(tmp)+
# 		geom_vline(xintercept = geos$time, color="lightgrey") +
# 	 	geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=median2, yend=median2), 
# 	 				 color="grey80", linewidth=2, lineend="round", alpha=0.1)+
# 	 	geom_ribbon(aes(x=timeMedian, ymin=HPD05median, ymax=HPD95median, color=NULL, group=grouping), fill=unique(tmp$colour), alpha=0.1)+
# 	 	geom_line(aes(x=timeMedian, y=median2, group=grouping, alpha=grouping), colour=unique(tmp$colour), linewidth=2, lineend="round")+
# 	 	scale_y_log10() + annotation_logticks(sides = 'l')+
# 	  facet_grid(data ~ root, scales="free")+
# 	 	scale_color_manual(values=unique(tmp$colour)) +
# 	 	scale_fill_manual(values=unique(tmp$colour)) +
# 	  scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
# 		scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 100),
# 						   minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
# 	 	labs(x="Time (Ma)", y="Median values",
# 	 		 title=paste0(i))+
# 	 	theme_classic()+
# 	 	theme(legend.position = "none")
# 	plot(fig3i)
# }
# dev.off()

#----
#---- Get extinction and mean diversification rates through time -----------------------------------

param$file = sub("_clads_.*", "", param$file)
rtt_raw$file = sub("_clads_.*", "", rtt_raw$file)

rates = merge(rtt_raw, param, by="file")
rates$div.x = NULL
rates$group.x = NULL
rates$calibration.x = NULL
rates$replicate.x = NULL
rates$root.x = NULL
colnames(rates) = c("file", "time", "rate", "HPD_05", "HPD_95", "mean", "revTime", "timePoint",
                   "calibration", "root", "replicate", "group", "div",
                   "sigma", "alpha", "epsilon", "lambda0", "meanLambdaTip", "logMeanLambdaTip", "m", "medianLambdaTip", "diversity")

# Extinction rate through time = epsilon*RTT (given than epsilon is extinction rate/speciation rate)
rates$extinction = rates$rate * rates$epsilon
rates$extinction05 = rates$HPD_05 * rates$epsilon
rates$extinction95 = rates$HPD_95 * rates$epsilon

# Net diversification rate = RTT*(1-epsilon)
rates$diversification = rates$rate * (1-rates$epsilon)
rates$diversification05 = rates$HPD_05 * (1-rates$epsilon)
rates$diversification95 = rates$HPD_95 * (1-rates$epsilon)

#
#---- Summarise all slopes -------------------------------------------------------------------------

ratess = rates %>% group_by(group, calibration, root, diversity, timePoint) %>% 
  summarise(timeMedian = median(revTime),
            timeHPD05 = hdi(revTime)[1],
            timeHPD95 = hdi(revTime)[2],
            speciation = median(rate),
            speciation05 = median(HPD_05),
            speciation95 = median(HPD_95),
            extinction = median(extinction),
            extinction05 = median(extinction05),
            extinction95 = median(extinction95),
            diversification = median(diversification),
            diversification05 = median(diversification05),
            diversification95 = median(diversification95))

#
#---- Prepare for plotting -------------------------------------------------------------------------

geos = subset(geo, time>min(ratess$timeHPD05))

{ratess$colour = ratess$group
  ratess$colour[which(ratess$colour=="Amoebozoa")]="royalblue1"
  ratess$colour[which(ratess$colour=="Nucletmycea")]="steelblue2"
  ratess$colour[which(ratess$colour=="Holozoa")]="steelblue4"
  ratess$colour[which(ratess$colour=="Metamonada")]="forestgreen"
  ratess$colour[which(ratess$colour=="Discoba")]="orange2"
  ratess$colour[which(ratess$colour=="Haptista")]="yellow1"
  ratess$colour[which(ratess$colour=="Cryptista")]="hotpink1"
  ratess$colour[which(ratess$colour=="Archaeplastida")]="darkseagreen3"
  ratess$colour[which(ratess$colour=="Rhizaria")]="darkorchid2"
  ratess$colour[which(ratess$colour=="Stramenopila")]="darkorchid3"
  ratess$colour[which(ratess$colour=="Alveolata")]="darkorchid4"}

ratess$group = factor(ratess$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

ratess$colour = factor(ratess$colour,
                    levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))

ratess$calibration = factor(ratess$calibration, levels=c("MC01", "MC02"))
ratess$root = factor(ratess$root, levels=c("rootD", "rootA"))

ratess$grouping = paste0(ratess$group, "-", ratess$diversity)
unique(ratess$grouping)
ratess$grouping = factor(ratess$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
                                                   "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
                                                   "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
                                                   "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))

write.table(ratess, paste0(files$outDir, "/", files$outPreffix, "rates.tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

#
#---- Plotting Speciation, extinction and net diversification rates --------------------------------

library(tidyr)
data = pivot_longer(ratess, -c(group, colour, calibration, root, diversity, grouping, timePoint, timeMedian, timeHPD05, timeHPD95))
data$rate = sub("05|95", "", data$name)
data$stat = ifelse(grepl("05", data$name), "HPD05", ifelse(grepl("95", data$name), "HPD95", "median"))
data$name = NULL
data = pivot_wider(data, names_from = "stat", values_from = "value")

data$rate = factor(data$rate, levels=c("diversification", "speciation", "extinction"))

(fig3 = ggplot(data)+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_segment(aes(x=timeHPD05, xend=timeHPD95, y=timeMedian, yend=timeMedian),
                 color="grey80", linewidth=2, lineend="round", alpha=0.1)+
    geom_ribbon(aes(x=timeMedian, ymin=HPD05, ymax=HPD95, fill=colour, color=NULL, group=grouping), alpha=0.1)+
    geom_line(aes(x=timeMedian, y=median, colour=group, alpha=grouping))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(rate ~ root + calibration, scales="free")+
    # facet_wrap(rate ~ root+diversity, scales="free_y", axes="all_y", 
    #            ncol=length(unique(c(data$diversity, data$root))), nrow=length(unique(data$rate)))+
    scale_color_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_fill_manual(values=as.character(sort(unique(rtt$colour)))) +
    scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
    scale_x_continuous(breaks=seq(min(round(rtt$timeMedian, -2)), 0, 200),
                       minor_breaks=seq(min(round(rtt$timeMedian, -2)), 0, 50)) +
    labs(x="Time (Ma)", y="Rates Through Time")+
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$outDir, "/", files$outPreffix, "rates.pdf"), width=11.69, height=8.27, paper='special'); plot(fig3); dev.off()


#----
# #---- And now plot all independent run -------------------------------------------------------------
# 
# rtt_raw$run = sub("\\/clads.*", "", rtt_raw$file)
# 
# {rtt_raw$colour = rtt_raw$group
# 	rtt_raw$colour[which(rtt_raw$colour=="Amoebozoa")]="royalblue1"
# 	rtt_raw$colour[which(rtt_raw$colour=="Nucletmycea")]="steelblue2"
# 	rtt_raw$colour[which(rtt_raw$colour=="Holozoa")]="steelblue4"
# 	rtt_raw$colour[which(rtt_raw$colour=="Metamonada")]="forestgreen"
# 	rtt_raw$colour[which(rtt_raw$colour=="Discoba")]="orange2"
# 	rtt_raw$colour[which(rtt_raw$colour=="Haptista")]="yellow1"
# 	rtt_raw$colour[which(rtt_raw$colour=="Cryptista")]="hotpink1"
# 	rtt_raw$colour[which(rtt_raw$colour=="Archaeplastida")]="darkseagreen3"
# 	rtt_raw$colour[which(rtt_raw$colour=="Rhizaria")]="darkorchid2"
# 	rtt_raw$colour[which(rtt_raw$colour=="Stramenopila")]="darkorchid3"
# 	rtt_raw$colour[which(rtt_raw$colour=="Alveolata")]="darkorchid4"}
# 
# rtt_raw$group = factor(rtt_raw$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
# 
# rtt_raw$colour = factor(rtt_raw$colour,
# 					 levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))
# 
# rtt_raw$grouping = paste0(rtt_raw$group, "-", rtt_raw$div)
# unique(ratess$grouping)
# rtt_raw$grouping = factor(rtt_raw$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
#                                                      "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
#                                                      "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
#                                                      "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))
# 
# pdf(paste0(files$outDir, "clads_RTT_allRuns.pdf"), width=11.69, height=8.27, paper='special')
# for(r in unique(rtt_raw$run)){
# 	cat("\r  ", grep(r, unique(rtt_raw$run)), "/", length(unique(rtt_raw$run)), "          ", sep="")
# 	tmp = subset(rtt_raw, run == r)
# 	geos = subset(geo, time>min(rtt_raw$revTime))
# 	plot = ggplot(tmp)+
# 	  geom_vline(xintercept = geos$time, color="lightgrey") +
# 	  geom_ribbon(aes(x=revTime, ymin=HPD_05, ymax=HPD_95, fill=colour, color=NULL, group=grouping), alpha=0.1)+
# 	  geom_line(aes(x=revTime, y=rate, colour=group, group=grouping, alpha=grouping))+
# 	  scale_y_log10() + annotation_logticks(sides = 'l')+
# 	  # facet_grid( ~ div, scales="free")+
# 	  scale_color_manual(values=as.character(sort(unique(tmp$colour)))) +
# 	  scale_fill_manual(values=as.character(sort(unique(tmp$colour)))) +
# 	  scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
# 	  scale_x_continuous(breaks=seq(min(round(tmp$revTime, -2)), 0, 100),
# 	                     minor_breaks=seq(min(round(tmp$revTime, -2)), 0, 50)) +
# 	  labs(title=paste0(r),
# 	       y="Mean speciation rate (Ln)", x="Time (Ma)")+
# 	  theme_classic()+
# 	  theme(legend.position = "none")
# 	plot(plot)
# }; rm(r, tmp, plot); dev.off()
# 
# # And the same for the DTT
# 
# dtt_raw$run = sub("\\/clads.*", "", dtt_raw$file)
# 
# {dtt_raw$colour = dtt_raw$group
# 	dtt_raw$colour[which(dtt_raw$colour=="Amoebozoa")]="royalblue1"
# 	dtt_raw$colour[which(dtt_raw$colour=="Nucletmycea")]="steelblue2"
# 	dtt_raw$colour[which(dtt_raw$colour=="Holozoa")]="steelblue4"
# 	dtt_raw$colour[which(dtt_raw$colour=="Metamonada")]="forestgreen"
# 	dtt_raw$colour[which(dtt_raw$colour=="Discoba")]="orange2"
# 	dtt_raw$colour[which(dtt_raw$colour=="Haptista")]="yellow1"
# 	dtt_raw$colour[which(dtt_raw$colour=="Cryptista")]="hotpink1"
# 	dtt_raw$colour[which(dtt_raw$colour=="Archaeplastida")]="darkseagreen3"
# 	dtt_raw$colour[which(dtt_raw$colour=="Rhizaria")]="darkorchid2"
# 	dtt_raw$colour[which(dtt_raw$colour=="Stramenopila")]="darkorchid3"
# 	dtt_raw$colour[which(dtt_raw$colour=="Alveolata")]="darkorchid4"}
# 
# dtt_raw$group = factor(dtt_raw$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
# 
# dtt_raw$colour = factor(dtt_raw$colour,
# 						 levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))
# 
# dtt_raw$grouping = paste0(dtt_raw$group, "-", dtt_raw$div)
# unique(ratess$grouping)
# dtt_raw$grouping = factor(dtt_raw$grouping, levels=c("Amoebozoa-max", "Amoebozoa-min", "Nucletmycea-max", "Nucletmycea-min", "Holozoa-max", "Holozoa-min",
#                                                      "Metamonada-max", "Metamonada-min", "Discoba-max", "Discoba-min",
#                                                      "Haptista-max", "Haptista-min", "Cryptista-max", "Cryptista-min", "Archaeplastida-max", "Archaeplastida-min",
#                                                      "Rhizaria-max", "Rhizaria-min", "Stramenopila-max", "Stramenopila-min", "Alveolata-max", "Alveolata-min"))
# 
# pdf(paste0(files$outDir, "clads_dtt_allRuns.pdf"), width=11.69, height=8.27, paper='special')
# for(r in unique(dtt_raw$run)){
# 	cat("\r  ", grep(r, unique(dtt_raw$run)), "/", length(unique(dtt_raw$run)), "          ", sep="")
# 	tmp = subset(dtt_raw, run == r)
# 	geos = subset(geo, time>min(dtt_raw$time))
# 	plot = ggplot(tmp)+
# 		geom_vline(xintercept = geos$time, color="lightgrey") +
# 		geom_ribbon(aes(x=time, ymin=min, ymax=max, fill=colour, color=NULL, group=grouping), alpha=0.1)+
# 		geom_line(aes(x=time, y=lineages, colour=group, group=grouping, alpha=grouping))+
# 		scale_y_log10() + annotation_logticks(sides = 'l')+
# 		# facet_grid( ~ div, scales="free")+
# 		scale_color_manual(values=as.character(sort(unique(tmp$colour)))) +
# 		scale_fill_manual(values=as.character(sort(unique(tmp$colour)))) +
# 	  scale_alpha_manual(values=c(0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1, 0.6, 1))+
# 		scale_x_continuous(breaks=seq(min(round(tmp$time, -2)), 0, 100),
# 						   minor_breaks=seq(min(round(tmp$time, -2)), 0, 50)) +
# 		labs(title=paste0(r),
# 			 y="Estimated diversity through time (Ln)", x="Time (Ma)")+
# 		theme_classic()+
# 		theme(legend.position = "none")
# 	plot(plot)
# }; rm(r, tmp, plot); dev.off()
# 
# 

#----
#----
# 
