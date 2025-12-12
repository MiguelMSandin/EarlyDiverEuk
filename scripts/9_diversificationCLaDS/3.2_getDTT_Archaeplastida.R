#----
#---- loading packages and data.frames -------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ape)
library(phangorn)
library(treeio)
library(tidytree)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----
#---- set working directory and file names ---------------------------------------------------------
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01-2/")

rm(list=ls()[!ls() %in% c("geo", "files")])
# .rs.restartR()

files = list(dtts=grep("_clads_DTT\\.tsv", dir(recursive=TRUE), value=TRUE),
             ltts="plots/LTT/MC01_LTT_Archaeplastida_summary.tsv",
             rdatas=grep("Archaeplastida.*\\.Rdata", dir(recursive=TRUE), value=TRUE),
             out="plots/clads/DTT_MC01-2_Archaeplastida") 
files

#----
#---- Extract DTT and subclades --------------------------------------------------------------------

# # just to understand what we are doing
for(i in 1:10){
 tree = rtree(10)
 plot(tree)
 cat("Choose a random node", i, "\n")
 treet = extract.clade(tree, interactive=TRUE)
 plot(treet)
 Sys.sleep(4)
};rm(i,tree,treet);cat("Done\n")

dtt = data.frame()
for(file in files$rdatas){
 # Open files
 load(file)
 root = sub("\\/.*", "", file)
 run = file %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)
 div = as.numeric(file %>% sub("_clads\\.Rdata", "", .) %>% sub(".*_div", "", .))
 cat("\n", root, " ", run, " ", div, " (", grep(file, files$rdatas), "/", length(files$rdatas), ")", sep="", "\n")
 tmp = data.frame()
 # Now loop through all replicates to extract (manually...) all subclades
 for(i in 1:10){
   tree = CladsOutput$enhanced_trees[[i]]$tree
   plot(tree)
   cat("  Select Rhodophyta...", i, " ")
   treet = extract.clade(tree, interactive=TRUE)
   tmpi=as.data.frame(ltt.plot.coords(treet))
   time_points = seq(min(tmpi$time),0, length.out=51)
   for(t in 1:51){
     ss = subset(tmpi, time<=time_points[t])
     tmp = rbind(tmp, data.frame(root=root, run=run, div=div, replicate=i, subclade="Rhodophyta",
                                 time_point=t, time=time_points[t], N=max(ss$N)))
   }
   cat("\r  Rhodophyta", i, "selected\n")
   cat("  Select Chloroplastida...", i, " ")
   treet = extract.clade(tree, interactive=TRUE)
   tmpi=as.data.frame(ltt.plot.coords(treet))
   time_points = seq(min(tmpi$time),0, length.out=51)
   for(t in 1:51){
     ss = subset(tmpi, time<=time_points[t])
     tmp = rbind(tmp, data.frame(root=root, run=run, div=div, replicate=i, subclade="Chloroplastida",
                                 time_point=t, time=time_points[t], N=max(ss$N)))
   }
   cat("\r  Chloroplastida", i, "selected\n")
 }
 dtt = rbind(dtt, tmp)
}
rm(file, root, run, div, tmp, i, tree, treet, tmpi, time_points, t, ss)

dtt$diversity = NA
unique(dtt$div)
tmp = mean(as.numeric(unique(dtt$div)))
dtt$diversity = fifelse(as.numeric(dtt$div) > tmp, "max", "min")

# Safe this hell of a table so we don't have to create it again
write.table(dtt, paste0(files$out, ".tsv"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
dtt = fread(paste0(files$out, ".tsv"))

#
#----
#----
#---- Read tables ----------------------------------------------------------------------------------

# Open extracted DTTs from Chloroplastida and Rhodophyta
tmp = fread(paste0(files$out, ".tsv"))
dtt = tmp %>% group_by(subclade, root, run, diversity, time_point) %>%
  summarise(time = median(time),
            lineages=median(N),
            min=min(N),
            max=max(N))
dtt$lineages = dtt$lineages + 1
dtt$min = dtt$min + 1
dtt$max = dtt$max + 1

# Open DTTs
tmp = data.frame()
for(file in files$dtts){
  if(grepl("(Discoba)|(Amoebozoa)|(Rhizaria)|(Alveolata)|(Nucletmycea)|(Holozoa)|(Archaeplastida)", file)){
    tmpi = fread(file)
    tmpi$time_points = 1:51
    tmpi$root = sub("\\/.*", "", file)
    tmpi$run = file %>% sub("root.\\/", "", .) %>% sub("\\/.*", "", .)
    tmpi$div = as.numeric(file %>% sub("_clads_DTT\\.tsv", "", .) %>% sub(".*_div", "", .))
    tmpi$clade = file %>% sub("\\/clade_.*", "", .) %>% sub(".*\\/", "", .)
    tmp = rbind(tmp, tmpi)
  }
};rm(tmpi)
tmp$diversity = NA
for(cl in unique(tmp$clade)){
  dm = mean(as.numeric(unique(subset(tmp, clade==cl)$div)))
  for(di in unique(subset(tmp, clade==cl)$div)){
      d = fifelse(di > dm, "max", "min")
      tmp$diversity[which(tmp$clade==cl & tmp$div==di)]=d
  }
};rm(cl, dm, di, d)

tmp = select(tmp, c(clade, root, run, diversity, time_points, time, lineages, min, max))
colnames(tmp) = c("clade", "root", "run", "diversity", "time_point", "time", "lineages", "min", "max")
colnames(dtt) = c("clade", "root", "run", "diversity", "time_point", "time", "lineages", "min", "max")

# bind datasets
all(names(dtt) == names(tmp))

dtt = rbind(dtt, tmp)

#
#----
#---- Summarise the table and prepapre for plotting ------------------------------------------------

summ = dtt %>% group_by(clade, root, diversity, time_point) %>% 
  summarise(t05=quantile(time, 0.05),
            t95=quantile(time, 0.95),
            time=median(time),
            n05=quantile(lineages, 0.05),
            n95=quantile(lineages, 0.95),
            n=median(lineages))

# Open LTTs
tmp = fread(files$ltts)
tmp = subset(tmp, replicates > max(tmp$replicates)*0.5)
tmp$n05 = tmp$N
tmp$n95 = tmp$N
tmp$diversity = "LTT"
tmp$time_points = NA
tmp = select(tmp, c(tree, root, diversity, time_points, q05, q95, median, n05, n95, N))
colnames(tmp) = c("clade", "root", "diversity", "time_point", "t05", "t95", "time", "n05", "n95", "n")

# bind 
summ = rbind(summ, tmp)

# Prettify
sort(unique(summ$clade))
summ$clade = factor(summ$clade, 
				   levels=c("Main_tree", "Archaeplastida", "Rhodophyta", "Chloroplastida",
				   		 "Discoba", "Rhizaria", "Alveolata", "Amoebozoa", "Nucletmycea", "Holozoa"))

{summ$colour = as.character(summ$clade)
	summ$colour[which(summ$colour=="Main_tree")]="black"
	summ$colour[which(summ$colour=="Archaeplastida")]="darkseagreen3"
	summ$colour[which(summ$colour=="Rhodophyta")]="red4"
	summ$colour[which(summ$colour=="Chloroplastida")]="olivedrab"
	summ$colour[which(summ$colour=="Chlorophyta")]="olivedrab3"
	summ$colour[which(summ$colour=="Streptophyta")]="olivedrab1"
	summ$colour[which(summ$colour=="Discoba")]="orange2"
	summ$colour[which(summ$colour=="Rhizaria")]="darkorchid2"
	summ$colour[which(summ$colour=="Alveolata")]="darkorchid4"
	summ$colour[which(summ$colour=="Amoebozoa")]="royalblue1"
	summ$colour[which(summ$colour=="Nucletmycea")]="steelblue2"
	summ$colour[which(summ$colour=="Holozoa")]="steelblue4"}

summ$colour = factor(summ$colour, 
				   levels=c("black", "darkseagreen3", "red4", "olivedrab", "olivedrab3", "olivedrab1",
				   		 "orange2", "darkorchid2", "darkorchid4", "royalblue1", "steelblue2", "steelblue4"))

summ$root = factor(summ$root, levels=c("rootD", "rootA"))
summ$diversity = factor(summ$diversity, levels=c("LTT", "max", "min"))

#
#----
#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table

geos = subset(geo, time > min(summ$t05, na.rm=TRUE))

(lttplot = ggplot(summ)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
    # geom_ribbon(aes(ymin=n05, ymax=n95, x=time, xmin=t05, xmax=t95, y=n, fill=clade), alpha=0.2, colour = NA) +
    geom_ribbon(aes(ymin=n05, ymax=n95, x=time, fill=clade), alpha=0.2, colour = NA) +
    geom_ribbon(aes(xmin=t05, xmax=t95, y=n, fill=clade), alpha=0.2, colour = NA) +
		geom_line(aes(x=time, y=n, color=clade)) +
		# geom_text(data=geos, aes(x=mid, y=max(summ$n)+5, label=period), size=2)+
		scale_x_continuous(breaks=seq((round(min(summ$t05, na.rm=TRUE), -2)), 0, 100),
						   minor_breaks=seq((round(min(summ$t05, na.rm=TRUE), -2)), 0, 50)) +
		scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(diversity ~ root, scales="free")+
		scale_color_manual(values=as.character(sort(unique(summ$colour))))+
		scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
		theme_classic()+
		labs(y="Ln (N)", x="Time (Ma)")+
    theme(legend.position = "none"))

pdf(paste0(files$out, ".pdf"), width=8.27, height=11.69, paper='special'); plot(lttplot); dev.off()



(lttplot = ggplot(subset(summ, clade=="Rhodophyta" | clade=="Chloroplastida" | clade=="Archaeplastida"))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    # geom_ribbon(aes(ymin=n05, ymax=n95, x=time, xmin=t05, xmax=t95, y=n, fill=clade), alpha=0.2, colour = NA) +
    geom_ribbon(aes(ymin=n05, ymax=n95, x=time, fill=clade), alpha=0.2, colour = NA) +
    # geom_ribbon(aes(xmin=t05, xmax=t95, y=n, fill=clade), alpha=0.2, colour = NA) +
    geom_line(aes(x=time, y=n, color=clade)) +
    # geom_text(data=geos, aes(x=mid, y=max(summ$n)+5, label=period), size=2)+
    scale_x_continuous(breaks=seq((round(min(summ$t05, na.rm=TRUE), -2)), 0, 100),
                       minor_breaks=seq((round(min(summ$t05, na.rm=TRUE), -2)), 0, 50)) +
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(diversity ~ root, scales="free")+
    scale_color_manual(values=as.character(sort(unique(summ$colour))))+
    scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
    theme_classic()+
    labs(y="Ln (N)", x="Time (Ma)")+
    theme(legend.position = "none"))


#----
