#---- 
#---- Loading packages -----------------------------------------------------------------------------

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

setwd("/home/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootA/")

files = list(summary=grep("medusa.*summary\\.tsv", dir(recursive=TRUE), value=TRUE),
			 shifts=grep("medusa.*shifts\\.tsv", dir(recursive=TRUE), value=TRUE),
			 outPrefix=paste0("../../../plots/medusa/", sub(".*\\/", "", getwd()), "_medusa"))

files

#---- 
#---- Reading files --------------------------------------------------------------------------------

files$clades = unique(files$shifts %>% sub(".*medusa\\/", "", .) %>% sub("_.*", "", .))
files$diversity = data.frame()
for(clade in files$clades){
	f = grep(clade, files$summary, value=TRUE)
	div = f %>% sub("_medusa.*", "", .) %>% sub(".*_d", "", .)
	div = as.numeric(div)
	m = mean(div)
	d = data.frame(file=f, diversity=fifelse(div > m, "max", "min"))
	files$diversity = rbind(files$diversity, d)
};rm(clade, f, div, m, d)

# Now read all files
summaries = data.frame()
for(file in files$summary){
	cat("\r  Loading file (", grep(file, files$summary), "/", length(files$summary), ") ", file, "                    " , sep="", end="")
	tmp = fread(file)
	div = files$diversity[grep(file, files$diversity[,1]),2]
	tmp$file = file
	tmp$group = file %>% sub(".*medusa\\/", "", .) %>% sub("_d.*", "", .)
	tmp$div = files$diversity[grep(file, files$diversity[,1]),2]
	summaries = rbind(summaries, tmp)
}; rm(tmp, div); cat("\n")

table(files$summary %>% sub(".*medusa\\/", "", .) %>% sub("_d.*", "", .))

#---- 
#---- Prepare for plotting -------------------------------------------------------------------------

# {summaries$colour = summaries$group
# 	summaries$colour[which(summaries$colour=="Amoebozoa")]="royalblue1"
# 	summaries$colour[which(summaries$colour=="Nucletmycea")]="steelblue2"
# 	summaries$colour[which(summaries$colour=="Holozoa")]="steelblue4"
# 	summaries$colour[which(summaries$colour=="Metamonada")]="forestgreen"
# 	summaries$colour[which(summaries$colour=="Discoba")]="orange2"
# 	summaries$colour[which(summaries$colour=="Haptista")]="yellow1"
# 	summaries$colour[which(summaries$colour=="Cryptista")]="hotpink1"
# 	summaries$colour[which(summaries$colour=="Archaeplastida")]="darkseagreen3"
# 	summaries$colour[which(summaries$colour=="Rhizaria")]="darkorchid2"
# 	summaries$colour[which(summaries$colour=="Stramenopila")]="darkorchid3"
# 	summaries$colour[which(summaries$colour=="Alveolata")]="darkorchid4"}
# 
summaries$group = factor(summaries$group, levels=c("Amoebozoa", "Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
# 
# summaries$colour = factor(summaries$colour,
# 					levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2", "darkorchid3","darkorchid4"))

#---- And plot -------------------------------------------------------------------------------------

geos = subset(geo, time>-max(summaries$age))

datas = subset(summaries, pass)
# datas = subset(summaries, count>2)
# datas = summaries
(shifts = ggplot(datas, aes(x=-age, y=group, fill=root))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_violin(scale="width", alpha=0.6, drop=FALSE)+
		# geom_point(aes(colour=root), alpha=0.4)+
		facet_wrap(~div, ncol=1)+
		scale_x_continuous(breaks=seq(-max(round(datas$age, -2)), 0, 100),
						   minor_breaks=seq(-max(round(datas$age, -2)), 0, 50)) +
		scale_fill_manual(values=c("turquoise4", "grey20")) +
		# scale_colour_manual(values=c("turquoise4", "grey20")) +
		labs(title="Shifts in diversification rates",
			 y="Clades", x="Time before present (Ma)")+
		theme_classic())
# geom_violin; scale: if "area" (default), all violins have the same area (before trimming the tails). If "count", areas are scaled proportionally to the number of observations. If "width", all violins have the same maximum width.
pdf(paste0(files$outPrefix, "_summaryShifts.pdf"), width=11.69, height=8.27, paper='special'); plot(shifts); dev.off()


# datas = subset(datas, age>540)
# geos = subset(geo, time < -min(datas$age) & time > -max(datas$age))
# (shiftsProt = ggplot(datas, aes(x=-age, y=group, fill=root))+
# 		geom_vline(xintercept = geos$time, color="lightgrey") +
# 		geom_violin(scale="width", alpha=0.6, drop=FALSE)+
# 		# geom_point(aes(colour=root), alpha=0.4)+
# 		facet_wrap(~div, ncol=1)+
# 		scale_x_continuous(breaks=seq(-max(round(datas$age, -2)), 0, 100),
# 						   minor_breaks=seq(-max(round(datas$age, -2)), 0, 50)) +
# 		scale_fill_manual(values=c("turquoise4", "grey20")) +
# 		# scale_colour_manual(values=c("turquoise4", "grey20")) +
# 		labs(title="Shifts in diversification rates",
# 			 y="Clades", x="Time before present (Ma)")+
# 		theme_classic())
# 
# pdf(paste0(files$outPrefix, "_summaryShifts_proterozoic.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsProt); dev.off()


#---- 
#---- Get how many shifts there are per clade ------------------------------------------------------

datas = summaries
datas = datas %>% group_by(group, div, file) %>% summarise(passed=sum(pass), total=length(unique(Shift.Node)))
datas = melt(as.data.table(datas), id.vars=c("group", "div", "file"))

ggplot(datas, aes(x=group, y=value))+
	geom_boxplot()+
	facet_wrap(variable ~ div, ncol=2, scales="free")+
	theme_classic()+
	theme(axis.text.x=element_text(angle=20, hjust=1))

ggplot(summaries, aes(x=count))+
	geom_histogram()+
	facet_wrap(group~div, scales="free_y")+
	theme_classic()+
	theme(axis.text.x=element_text(angle=20, hjust=1))

datas = summaries %>% group_by(group, div, file) %>% summarise(c1=sum(count == 1), c2=sum(count <= 2), c3=sum(count <= 3), c4=sum(count <= 4), c5=sum(count <= 5),
															   c6=sum(count <= 6), c7=sum(count <= 7), c8=sum(count <= 8), c9=sum(count <= 9), c10=sum(count <= 10))
datas = melt(as.data.table(datas), id.vars=c("group", "div", "file"))
datas = datas %>% group_by(group, div, variable) %>% summarise(min=min(value), q05=quantile(value, 0.05), q25=quantile(value, 0.25),
																   mean=mean(value),
																   q75=quantile(value, 0.75), q95=quantile(value, 0.95), max=max(value))

ggplot(datas, aes(x=variable, color=group))+
	# geom_ribbon(aes(ymin=q05, ymax=q95, fill=group), alpha=0.4) +
	geom_line(aes(y=q05, group=group), alpha=0.2)+
	geom_line(aes(y=q95, group=group), alpha=0.2)+
	geom_line(aes(y=mean, group=group))+
	facet_wrap(~div)+
	theme_classic()

#---- 
#---- Others ---------------------------------------------------------------------------------------

# div = 0.4
# reps = 10000
# tmp = data.table(x=1:reps, n=rnorm(reps, 1/div, 0.5), r=rep(0,reps), t=rep(0,reps))
# for(i in 1:reps){
# 	cat("\r  ", round(i/reps*100), "%", sep="")
# 	r=0
# 	while(r < 1){r = round(rnorm(1, 1/div, 0.5))}
# 	tmp$t[i] = r
# 	tmp$r[i] = round(tmp$n[i])
# };rm(i, r);cat("\n")
# 
# cat("  There are ", reps/div, " theoretical tips and ", sum(tmp$t)," randomly generated tips\n", sep="")
# 
# tmpm = melt(tmp, id.vars=c("x"))
# ggplot(tmpm, aes(value, fill=variable))+
# 	geom_histogram(position = "identity", alpha=0.6)+
# 	theme_minimal()

#---- 


