#---- 
#---- Loading packages -----------------------------------------------------------------------------

library(ape)
library(dplyr)
library(data.table)
library(coda)
library(HDInterval)
library(ggplot2)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				 era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#---- 
#---- Set working directory and file names ---------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo")])
# .rs.restartR()

setwd("/home/miguel/Desktop/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootA/")

files = list(dirFiles="treePar",
			 dirTrees="clades",
			 dirOut="../../../plots/treePar/",
			 # outPrefix="")
			 outPrefix=sub(".*\\/", "", getwd()))

files$rates = grep(paste0(files$dirFiles, ".*\\/.*\\/.*", "EBD_rates.tsv"), dir(recursive=TRUE), value=TRUE)
files$replicates = unique(sub("\\/.*", "", files$rates))
files$clades = unique(files$rates %>% sub(paste0(".*", files$dirFiles, "\\/"), "", .) %>% sub("\\/.*", "", .))

files

#---- 
#---- Read files -----------------------------------------------------------------------------------

rates = data.table(); i = 0; count = list()
for(cl in files$clades){
	count[cl] = 0
}
for(re in files$replicates){
	filesre = grep(re, files$rates, value=TRUE)
	for(cl in files$clades){
		if(any(grepl(cl, filesre))){
			count[[cl]] = count[[cl]] + 1
			tmp = grep(cl, filesre, value=TRUE)
			intervals = as.numeric(tmp %>% gsub(".*_i", "", .) %>% gsub("_.*", "", .))
			tmp = fread(tmp)
			tmp$replicate = re
			tmp$clade = cl
			tmp$time_unit = rep(0:intervals, 2)
			rates = rbind(rates, tmp)
		}
	}
}; rm(re, cl, filesre, tmp, intervals)

count

#---- 
#---- Transform and summarise tables for plotting --------------------------------------------------

tmpr = select(rates, c("div", "time", "time_unit", "speciation", "extinction", "diversification", "diversity", "clade", "replicate"))
tmphi = select(rates, c("div", "time", "time_unit", "spec_HPDmin", "ext_HPDmin", "div_HPDmin", "diversity", "clade", "replicate"))
tmpha = select(rates, c("div", "time", "time_unit", "spec_HPDmax", "ext_HPDmax", "div_HPDmax", "diversity", "clade", "replicate"))

tmpr = melt(tmpr, id.vars=c("replicate", "clade", "div", "diversity", "time", "time_unit"))
tmphi = melt(tmphi, id.vars=c("replicate", "clade", "div", "diversity", "time", "time_unit"))
tmpha = melt(tmpha, id.vars=c("replicate", "clade", "div", "diversity", "time", "time_unit"))

colnames(tmpr) = c("replicate", "clade", "div", "diversity", "time", "time_unit", "variable", "rate")
colnames(tmphi) = c("replicate", "clade", "div", "diversity", "time", "time_unit", "variable", "HPD_min")
colnames(tmpha) = c("replicate", "clade", "div", "diversity", "time", "time_unit", "variable", "HPD_max")

tmphi$variable = fifelse(grepl("spec", tmphi$variable), "speciation", 
						 fifelse(grepl("ext", tmphi$variable), "extinction", 
						 		fifelse(grepl("div", tmphi$variable), "diversification", NA)))
tmpha$variable = fifelse(grepl("spec", tmpha$variable), "speciation", 
						 fifelse(grepl("ext", tmpha$variable), "extinction", 
						 		fifelse(grepl("div", tmpha$variable), "diversification", NA)))

all(tmpr$div == tmphi$div & tmpr$div == tmpha$div)
all(tmpr$diversity == tmphi$diversity & tmpr$diversity == tmpha$diversity)
all(tmpr$time == tmphi$time & tmpr$time == tmpha$time)
all(tmpr$time_unit == tmphi$time_unit & tmpr$time_unit == tmpha$time_unit)
all(tmpr$clade == tmphi$clade & tmpr$clade == tmpha$clade)
all(tmpr$variable == tmphi$variable & tmpr$variable == tmpha$variable)

ratesl = tmpr
ratesl$HPDmin = tmphi$HPD_min
ratesl$HPDmax = tmpha$HPD_max
rm(tmpr, tmphi, tmpha)

# Summarize ________________________________________________________________________________________

ratess = ratesl %>% group_by(clade, diversity, variable, time_unit) %>% summarise(timeHPDmin=hdi(time)[1],
																				  timeHPDmax=hdi(time)[2],
																				  time=median(time),
																				  # rateMedianHPDmin=hdi(rate, credMass=0.90)[1],
																				  # rateMedianHPDmax=hdi(rate, credMass=0.90)[2],
																				  rateMedianQ25=quantile(rate, 0.25),
																				  rateMedianQ75=quantile(rate, 0.75),
																				  rate=median(rate),
																				  rateMinMedian=median(HPDmin), 
																				  rateMaxMedian=median(HPDmax))

# Now discretize the rates by time intervals _______________________________________________________
ratesd = data.table()
for(cl in unique(ratess$clade)){
	cat("\r  Clade ", cl, " (", grep(cl, unique(ratess$clade)), "/", length(unique(ratess$clade)), ")          ", sep="")
	ss = subset(ratess, clade==cl)
	for(d in unique(ss$diversity)){
		ss2 = subset(ss, diversity==d)
		inter = (ss2$time[2]-ss2$time[1])*0.999
		for(i in 1:nrow(ss2)){
			ss3 = ss2[i,]
			ratesd = rbind(ratesd, ss3)
			ss3$time = ss3$time + inter
			ratesd = rbind(ratesd, ss3)
		}
	}
}; rm(d, ss, ss2, ss3, inter, i)

{ratesd$colour = ratesd$clade
	ratesd$colour[which(ratesd$colour=="Amoebozoa")]="royalblue1"
	ratesd$colour[which(ratesd$colour=="Nucletmycea")]="steelblue2"
	ratesd$colour[which(ratesd$colour=="Holozoa")]="steelblue4"
	ratesd$colour[which(ratesd$colour=="Metamonada")]="forestgreen"
	ratesd$colour[which(ratesd$colour=="Discoba")]="orange2"
	ratesd$colour[which(ratesd$colour=="Haptista")]="yellow1"
	ratesd$colour[which(ratesd$colour=="Cryptista")]="hotpink1"
	ratesd$colour[which(ratesd$colour=="Archaeplastida")]="darkseagreen3"
	ratesd$colour[which(ratesd$colour=="Rhizaria")]="darkorchid2"
	ratesd$colour[which(ratesd$colour=="Stramenopila")]="darkorchid3"
	ratesd$colour[which(ratesd$colour=="Alveolata")]="darkorchid4"}

ratesd$variable = factor(ratesd$variable, levels=c("speciation", "extinction", "diversification"))
ratesd$clade = factor(ratesd$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

ratesd$colour = factor(ratesd$colour,
					levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2","darkorchid3","darkorchid4"))

ratesd$discrete = rep(c(0, 1), nrow(ratesd)/2)

#---- 
#---- Plot -----------------------------------------------------------------------------------------

geos = subset(geo, time>-max(ratesd$time))

# (ratesplot = ggplot(ratesd)+
(ratesplot = ggplot(subset(ratesd, discrete==0))+
# (ratesplot = ggplot(subset(ratesd, discrete==0 & variable=="diversification"))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_ribbon(aes(x=-time, ymin=rateMedianQ25, ymax=rateMedianQ75, color=NULL, fill=clade), alpha=0.1)+
		geom_line(aes(x=-time, y=rate, colour=clade))+
		# geom_smooth(aes(x=-time, y=rate, colour=clade), se=FALSE)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq(-max(round(ratesd$time, -2)), 0, 200),
						   minor_breaks=seq(-max(round(ratesd$time, -2)), 0, 50)) +
		labs(y="Mean rates (Ln)", x="Time (Ma)", title=paste0(files$outPrefix))+
		facet_grid(variable~diversity, scales="free")+
		scale_color_manual(values=as.character(sort(unique(ratesd$colour)))) +
		scale_fill_manual(values=as.character(sort(unique(ratesd$colour)))) +
		theme_classic()+
		theme(legend.position = "none"))

# And export the table and the plot
write.table(ratess, paste0(files$dirOut, files$outPrefix, "_EBD_rates.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
pdf(paste0(files$dirOut, files$outPrefix, "_EBD_rates.pdf"), width=8.27, height=11.69, paper='special'); plot(ratesplot); dev.off()

#---- 
#---- Plot individual slopes -----------------------------------------------------------------------

pdf(paste0(files$dirOut, files$outPrefix, "_EBD_rates_clades.pdf"), width=8.27, height=11.69, paper='special')
for(cl in unique(ratesl$clade)){
	data = subset(ratesl, clade==cl)
	datas = data %>% group_by(diversity, variable, time_unit) %>% summarise(rate_mean=mean(rate), rate_median=median(rate),
																			time=median(time))
	geos = subset(geo, time>-max(data$time))
	ratesPlotClades = ggplot(data)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_line(aes(x=-time, y=rate, colour=replicate))+
		geom_line(data=datas, aes(x=-time, y=rate_mean), colour="orange")+
		geom_line(data=datas, aes(x=-time, y=rate_median), colour="orangered4")+
		# geom_smooth(aes(x=-time, y=rate, colour=clade), se=FALSE)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq(-max(round(rates$time, -2)), 0, 200),
						   minor_breaks=seq(-max(round(rates$time, -2)), 0, 50)) +
		labs(y="Mean rates (Ln)", x="Time (Ma)", title=paste0(files$outPrefix, ": ", unique(data$clade)))+
		facet_grid(variable~diversity, scales="free")+
		scale_color_manual(values=rep("grey60", length(unique(data$replicate)))) +
		theme_classic()+
		theme(legend.position = "none")
	plot(ratesPlotClades)
}; dev.off()

#---- 
#---- Plot different clades together ---------------------------------------------------------------

data = subset(ratesd, (clade=="Holozoa" | clade=="Rhizaria") & discrete==0)
# data = subset(ratesd, (clade=="Archaeplastida" | clade=="Cryptista" | clade=="Stramenopila" | clade=="Alveolata" | clade=="Haptista") & discrete==0)
geos = subset(geo, time>-max(data$time))
means = ratesd %>% group_by(variable, diversity) %>% summarise(median=median(rate), mean=mean(rate))
(ggplot(data)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_hline(data=means, aes(yintercept = means$mean), color="orangered4") +
		geom_ribbon(aes(x=-time, ymin=rateMedianQ25, ymax=rateMedianQ75, color=NULL, fill=clade), alpha=0.1)+
		geom_line(aes(x=-time, y=rate, colour=clade))+
		# geom_smooth(aes(x=-time, y=rate, colour=clade), se=FALSE)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 200),
						   minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
		labs(y="Mean rates (Ln)", x="Time (Ma)", title=paste0(files$outPrefix))+
		facet_grid(variable~diversity, scales="free")+
		scale_color_manual(values=as.character(sort(unique(data$colour)))) +
		scale_fill_manual(values=as.character(sort(unique(data$colour)))) +
		theme_classic())

#---- 
