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

# setwd("/home/miguel/Desktop/miguel/Documents/PhD/0_Thesis/2_chapter/data/diver/ebd/")
# setwd("/home/miguel/Desktop/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/rRAng1-2-1/treePar/")
setwd("/home/miguel/Desktop/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/fRAcat2-2/treeParL/")


files = list(clade="Holozoa",
			 dirTrees="../clades/")

{files$clade = gsub("\\/+$", "", files$clade)
	files$speciation = grep("speciation", dir(files$clade), value=TRUE)
	files$extinction = grep("extinction", dir(files$clade), value=TRUE)
	files$intervals = as.numeric(unique(c(files$speciation, files$extinction) %>% gsub(".*_i", "", .) %>% gsub("_.*", "", .)))
	files$diversities = as.numeric(unique(c(files$speciation, files$extinction) %>% gsub(".*_d", "", .) %>% gsub("_.*", "", .)))
	files$tree = paste0(files$dirTrees, grep(paste0(files$clade, "\\.tre"), dir(files$dirTrees), value=TRUE))}

if(length(files$intervals) != 1){cat("  Warning! Different files have different time intervals...")}
if(length(files$diversities) != 2){cat("  Warning! There should be only 2 different diversity fractions...")}

files

#---- 
#---- Read files -----------------------------------------------------------------------------------

tree = read.tree(files$tree)

files$root_age = -min(ltt.plot.coords(tree))
files$times = c(0, files$root_age * (1:files$intervals) / (files$intervals))

chains = 100001
rateS = data.table()
for(i in files$speciation){
	tmp = fread(paste0(files$clade, "/", i))
	tmp$file = i
	chains = c(chains, nrow(tmp))
	rateS = rbind(rateS, tmp)
}
rateE = data.table()
for(i in files$extinction){
	tmp = fread(paste0(files$clade, "/", i))
	tmp$file = i
	chains = c(chains, nrow(tmp))
	rateE = rbind(rateE, tmp)
}

# chains = chains[-1]
if(var(chains) != 0){cat("\n  Warning! diversity files have different number of generations...\n\n")}else{chains=unique(chains)}

#---- 
#---- Check convergence ----------------------------------------------------------------------------

ess = data.frame()
for(i in unique(rateS$file)){
	tmp=subset(rateS, file==i)
	plot(tmp$Likelihood)
	title(i)
	tmp = apply(select(tmp, -file), 2, effectiveSize)
	tmp = as.data.frame(tmp)
	tmp$file = i
	tmp$variable = row.names(tmp)
	ess = rbind(ess, tmp)
}
for(i in unique(rateE$file)){
	tmp=subset(rateE, file==i)
	tmp = apply(select(tmp, -file), 2, effectiveSize)
	tmp = as.data.frame(tmp)
	tmp$file = i
	tmp$variable = row.names(tmp)
	ess = rbind(ess, tmp)
}

ess$time.interval = NA
for(i in 1:nrow(ess)){
	if(grepl("\\[", ess$variable[i])){
		ess$time.interval[i] = as.numeric(ess$variable[i] %>% sub(".*\\[", "", .) %>% sub("\\]", "", .))-1
	}
}; rm(i, tmp)

ggplot(subset(ess, is.na(time.interval)), aes(x=variable, y=tmp, colour=file))+
	geom_point()+
	geom_hline(yintercept = 200, color="orangered4")+
	geom_hline(yintercept = 100, color="orangered1")+
	theme_classic()+
	theme(legend.position = "bottom")

# ggplot(subset(ess, !is.na(time.interval)), aes(x=-time.interval, y=tmp, colour=file))+
# 	geom_point()+
# 	geom_hline(yintercept = 200, color="orangered4")+
# 	geom_hline(yintercept = 100, color="orangered1")+
# 	theme_classic()+
# 	theme(legend.position = "bottom")

ess$ess2 = ifelse(ess$tmp > 400, 400, ess$tmp)
ggplot(subset(ess, !is.na(time.interval)), aes(x=-time.interval, y=ess2, colour=file))+
	geom_point()+
	geom_hline(yintercept = 200, color="orangered1")+
	geom_hline(yintercept = 100, color="orangered3")+
	geom_hline(yintercept = 50, color="orangered4")+
	theme_classic()+
	theme(legend.position = "bottom")

#---- 
#---- Transform and summarise tables for plotting --------------------------------------------------

rateS = melt(rateS, id.vars=c("file", "Iteration", "Posterior", "Likelihood", "Prior"))
rateE = melt(rateE, id.vars=c("file", "Iteration", "Posterior", "Likelihood", "Prior"))

rateS$div = as.numeric(rateS$file %>% gsub(".*_d", "", .) %>% gsub("_.*", "", .))
rateE$div = as.numeric(rateE$file %>% gsub(".*_d", "", .) %>% gsub("_.*", "", .))

rateS$time = NA
for(i in 1:length(unique(rateS$variable))){
	cat("\r  ", round(i/length(unique(rateS$variable))*100), "%", sep="")
	rateS$time[which(rateS$variable==unique(rateS$variable)[i])] = files$times[i]
}
rateE$time = NA
for(i in 1:length(unique(rateE$variable))){
	cat("\r  ", round(i/length(unique(rateE$variable))*100), "%", sep="")
	rateE$time[which(rateE$variable==unique(rateE$variable)[i])] = files$times[i]
};rm(i)

all(rateE$Iteration == rateS$Iteration)
all(rateE$div == rateS$div)
rateD = data.table(div=rateS$div,
				   iteration=rateS$Iteration,
				   time=rateS$time,
				   speciation=rateS$value,
				   extinction=rateE$value,
				   diversification=rateS$value-rateE$value)

rateDs = rateD %>% group_by(div, time) %>% summarise(spec_HPDmin=hdi(speciation)[[1]], spec_HPDmax=hdi(speciation)[[2]],
													 speciation=median(speciation), 
													 ext_HPDmin=hdi(extinction)[[1]], ext_HPDmax=hdi(extinction)[[2]],
													 extinction=median(extinction),
													 div_HPDmin=hdi(diversification)[[1]], div_HPDmax=hdi(diversification)[[2]],
													 diversification=median(diversification))
rateDs$diversity = NA
for(i in 1:nrow(rateDs)){
	rateDs$diversity[i] = ifelse(rateDs$div[i] < mean(rateDs$div), "min", "max")
}
rateDs

# Now discretize the rates by time intervals
rates = data.table()
for(d in unique(rateDs$diversity)){
	ss = subset(rateDs, diversity==d)
	inter = (ss$time[2]-ss$time[1])*0.999
	for(i in 1:nrow(ss)){
		sss = ss[i,]
		rates = rbind(rates, sss)
		sss$time = sss$time + inter
		rates = rbind(rates, sss)
	}
};rm(d, ss, inter, i)


#---- 
#---- Plot -----------------------------------------------------------------------------------------

geos = subset(geo, time>-files$root_age)

(ratesplot = ggplot(rates)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		# geom_ribbon(aes(x=-time, ymin=spec_HPDmin, ymax=spec_HPDmax, color=NULL), fill="steelblue2", alpha=0.1)+
		# geom_ribbon(aes(x=-time, ymin=ext_HPDmin, ymax=ext_HPDmax, color=NULL), fill="orangered2", alpha=0.1)+
		# geom_ribbon(aes(x=-time, ymin=div_HPDmin, ymax=div_HPDmax, color=NULL), fill="turquoise4", alpha=0.1)+
		geom_line(aes(x=-time, y=speciation), colour="steelblue2")+
		geom_line(aes(x=-time, y=extinction), colour="orangered2")+
		geom_line(aes(x=-time, y=diversification), colour="turquoise4")+
		# geom_smooth(aes(x=-time, y=diversification), method="lm", se=FALSE)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq(-max(round(rates$time, -2)), 0, 100),
						   minor_breaks=seq(-max(round(rates$time, -2)), 0, 50)) +
		labs(title=paste0(files$clade, "_i", files$intervals, "_EBD_rates"),
			 y="Mean rates (Ln)", x="Time (Ma)")+
		facet_wrap(~diversity, ncol=1, scales="free")+
		theme(legend.position = "none")+
		theme_classic())

# pdf(paste0(files$clade, "/", files$clade, "_i", files$intervals, "_EBD_rates.pdf"), width=11.69, height=8.27, paper='special'); plot(ratesplot); dev.off()
# write.table(rateDs, paste0(files$clade, "/", files$clade, "_i", files$intervals, "_EBD_rates.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#---- 
