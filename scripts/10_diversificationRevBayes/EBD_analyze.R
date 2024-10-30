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

setwd("/home/miguel/Desktop/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/rRAgamma1-1/")

files = list(dirFiles="treePar/",
			 dirTrees="clades/")

files$speciation = grep("speciation", dir(files$dirFiles, recursive=TRUE), value=TRUE)
files$extinction = grep("extinction", dir(files$dirFiles, recursive=TRUE), value=TRUE)
files$clades = unique(sub("\\/.*", "", c(files$speciation, files$extinction)))

files

#---- 
#---- Read files -----------------------------------------------------------------------------------

rates = data.table()
for(cl in files$clade){
	file = grep("EBD_rates.tsv", dir(paste0(files$dirFiles, "/", cl)), value=TRUE)
	tmp = fread(paste0(files$dirFiles, "/", cl, "/", file))
	tmp$clade = cl
	rates = rbind(rates, tmp)
}; rm(cl, file, tmp)

#---- 
#---- Transform and summarise tables for plotting --------------------------------------------------

tmpr = select(rates, c("div", "time", "speciation", "extinction", "diversification", "diversity", "clade"))
tmphi = select(rates, c("div", "time", "spec_HPDmin", "ext_HPDmin", "div_HPDmin", "diversity", "clade"))
tmpha = select(rates, c("div", "time", "spec_HPDmax", "ext_HPDmax", "div_HPDmax", "diversity", "clade"))

tmpr = melt(tmpr, id.vars=c("clade", "div", "diversity", "time"))
tmphi = melt(tmphi, id.vars=c("clade", "div", "diversity", "time"))
tmpha = melt(tmpha, id.vars=c("clade", "div", "diversity", "time"))

colnames(tmpr) = c("clade", "div", "diversity", "time", "variable", "rate")
colnames(tmphi) = c("clade", "div", "diversity", "time", "variable", "HPD_min")
colnames(tmpha) = c("clade", "div", "diversity", "time", "variable", "HPD_max")

tmphi$variable = fifelse(grepl("spec", tmphi$variable), "speciation", 
						 fifelse(grepl("ext", tmphi$variable), "extinction", 
						 		fifelse(grepl("div", tmphi$variable), "diversification", NA)))
tmpha$variable = fifelse(grepl("spec", tmpha$variable), "speciation", 
						 fifelse(grepl("ext", tmpha$variable), "extinction", 
						 		fifelse(grepl("div", tmpha$variable), "diversification", NA)))

all(tmpr$div == tmphi$div & tmpr$div == tmpha$div)
all(tmpr$diversity == tmphi$diversity & tmpr$diversity == tmpha$diversity)
all(tmpr$time == tmphi$time & tmpr$time == tmpha$time)
all(tmpr$clade == tmphi$clade & tmpr$clade == tmpha$clade)
all(tmpr$variable == tmphi$variable & tmpr$variable == tmpha$variable)

ratesl = tmpr
ratesl$HPDmin = tmphi$HPD_min
ratesl$HPDmax = tmpha$HPD_max
rm(tmpr, tmphi, tmpha)

# Now discretize the rates by time intervals
ratesd = data.table()
for(cl in unique(ratesl$clade)){
	ss = subset(ratesl, clade==cl)
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
		geom_vline(xintercept = geos$time, color="lightgrey") +
		# geom_ribbon(aes(x=-time, ymin=HPDmin, ymax=HPDmax, color=NULL, fill=clade), alpha=0.1)+
		geom_line(aes(x=-time, y=rate, colour=clade), linewidth=1)+
		# geom_smooth(aes(x=-time, y=rate, colour=clade), se=FALSE)+
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_x_continuous(breaks=seq(-max(round(ratesd$time, -2)), 0, 100),
						   minor_breaks=seq(-max(round(ratesd$time, -2)), 0, 50)) +
		labs(y="Mean rates (Ln)", x="Time (Ma)")+
		facet_grid(diversity~variable, scales="free")+
		scale_color_manual(values=as.character(sort(unique(ratesd$colour)))) +
		scale_fill_manual(values=as.character(sort(unique(ratesd$colour)))) +
		theme(legend.position = "none")+
		theme_classic())

# And export the table and the plot
write.table(ratesd, paste0(files$dirFiles, "EBD_rates.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
pdf(paste0(files$dirFiles, "EBD_rates.pdf"), width=11.69*2, height=8.27, paper='special'); plot(ratesplot); dev.off()

#---- 
