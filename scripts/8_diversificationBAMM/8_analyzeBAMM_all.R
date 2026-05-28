#---- 
#---- loading packages and functions ---------------------------------------------------------------

library(data.table)
library(dplyr)
library(ape)
library(paleotree)
library(phangorn)
library(ggplot2)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
                 timeLow= c(-2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58, 0),
                 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                 era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"),
                 eon=c("Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Proterozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic","Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic", "Phanerozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)


getBLbin = function(tree, max=Inf, min=0, geoBin=NA, verbose=FALSE){
  if(max == Inf){
    max = max(node.depth.edgelength(tree))
  }
  if(!is.na(geoBin)){
    if(geoBin %in% geo$eon){
      max=max(-subset(geo, eon==geoBin)$time)
      min=min(-subset(geo, eon==geoBin)$timeLow)
    }else if(geoBin %in% geo$era){
      max=max(-subset(geo, era==geoBin)$time)
      min=min(-subset(geo, era==geoBin)$timeLow)
    }else if(geoBin %in% geo$period){
      max=max(-subset(geo, period==geoBin)$time)
      min=min(-subset(geo, period==geoBin)$timeLow)
    }else{
      cat("  Error! 'geoBin' not found in the geological scale names.\n")
      stop()
    } 
  }
  if(is.null(tree$root.time)){tree$root.time = max(node.depth.edgelength(tree))}
  if(min >= tree$root.time){
    out = NA
  }else if(max >= tree$root.time & min == 0){
    out = sum(tree$edge.length)
  }else if(max >= tree$root.time & min != 0){
    trees = timeSliceTree(tree, min, plot=FALSE)
    out = sum(trees$edge.length)
  }else if(max < tree$root.time & min == 0){
    trees = timeSliceTree(tree, max, plot=FALSE)
    out = sum(tree$edge.length) - sum(trees$edge.length)
  }else{
    trees = timeSliceTree(tree, min, plot=FALSE)
    out = sum(trees$edge.length)
    trees = timeSliceTree(tree, max, plot=FALSE)
    out = out - sum(trees$edge.length)
  }
  if(verbose){
    return(cat("Max date: ", max, " - Min date: ", min, " (time bin: ", max-min, "); total branch length:\n", out, sep="", end="\n"))
  }else{
    return(out)
  }
}


#---- 
# 'setwd' to the root directory where you have RTT and shift files
setwd("/home/miguel/Documents/Uppsala/1_ecoEvo/data/euk/stepDating/")
rm(list=ls()[!ls() %in% c("files", "geo", "getBLbin")])
# .rs.restartR()
#---- 
#---- Set file names -------------------------------------------------------------------------------

files = list(files=grep("bamm\\/.*RTT\\.tsv$", dir(recursive=TRUE), value=TRUE),
             shifts=grep("bamm\\/.*shifts\\.tsv$", dir(recursive=TRUE), value=TRUE), 
             trees=grep("clades\\/clade_[A-Z][a-z]+\\.tre$", dir(recursive=TRUE), value=TRUE),
             dirPlots="plots/shifts/",
             RTTplot=paste0("bamm_RTT_allMedian_", sub(".*\\/", "", getwd()), ".pdf"),
             shiftsPlot=paste0("bamm_shifts_", sub(".*\\/", "", getwd())),
             shiftsData=paste0("bamm_shifts_", sub(".*\\/", "", getwd()), "_rawData.tsv"),
             bl=paste0("bamm_shifts_", sub(".*\\/", "", getwd()), "_branchLength.tsv"))

files$trees = grep("_annot", files$trees, invert=TRUE, value=TRUE)

if(!dir.exists(files$dirPlots)){dir.create(files$dirPlots)}

# data = fread(paste0(files$dirPlots, "/", files$shiftsData))
# bl = fread(paste0(files$dirPlots, "/", files$bl))
# data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
# data$change = factor(data$change, levels=c("root", "growth", "decay"))
# data$root = factor(data$root, levels=c("rootD", "rootA"))

#---- 
# #---- plot RTT BAMM output -------------------------------------------------------------------------
# 
# rm(list=ls()[!ls() %in% c("geo", "files")])
# 
# # Reading all RTT files
# data = data.frame()
# for(file in files$files){
#   cat("\r  Working on file (", grep(file, files$files), "/", length(files$files), ") ", file, "          ", sep="")
#   tmp = fread(file)
#   tmp$file = file
#   tmp$root = sub("\\/.*", "", file)
#   tmp$clade = file %>% sub(".*diversification_", "", .) %>% sub("_.*", "", .)
#   tmp$diversity = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
#   data = rbind(data, tmp)
# };rm(file, tmp); cat("\n")
# 
# # Create another variable for minimum and maximum diversity
# data$div = NA
# for(clad in unique(data$clade)){
#   ss = subset(data, clade==clad)
#   m = mean(unique(ss$diversity))
#   data$div[which(data$clade==clad & data$diversity < m)] = "min"
#   data$div[which(data$clade==clad & data$diversity > m)] = "max"
# }; rm(clad, ss, m)
# 
# # Create a variable for the time point sampled per file
# data$timePoint = rep(seq(1:100), length(unique(data$file))*length(unique(data$rate)))
# 
# # Summarise median slope per clade
# summ = data %>% group_by(root, clade, div, rate, timePoint) %>% summarise(timeMedian=median(time),
#                                                                           timeSD=sd(time),
#                                                                           median2=median(median),
#                                                                           HPDlowMedian=median(HPDlow),
#                                                                           HPDuppMedian=median(HPDupp))
# 
# # Adding colours
# {summ$colour = summ$clade
#   summ$colour[which(summ$colour=="Amoebozoa")] =      "royalblue1"
#   summ$colour[which(summ$colour=="Nucletmycea")] =    "steelblue2"
#   summ$colour[which(summ$colour=="Holozoa")] =        "steelblue4"
#   summ$colour[which(summ$colour=="Metamonada")] =     "forestgreen"
#   summ$colour[which(summ$colour=="Discoba")] =        "orange2"
#   summ$colour[which(summ$colour=="Haptista")] =       "yellow1"
#   summ$colour[which(summ$colour=="Cryptista")] =      "hotpink1"
#   summ$colour[which(summ$colour=="Archaeplastida")] = "darkseagreen3"
#   summ$colour[which(summ$colour=="Rhizaria")] =       "darkorchid2"
#   summ$colour[which(summ$colour=="Stramenopila")] =   "darkorchid3"
#   summ$colour[which(summ$colour=="Alveolata")] =      "darkorchid4"}
# 
# # Adding factors for plotting
# summ$rate = factor(summ$rate, levels=c("extinction","speciation","diversification"))
# summ$root = factor(summ$root, levels=c("rootD","rootA"))
# summ$clade = factor(summ$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
# summ$colour = factor(summ$colour, levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2","darkorchid3","darkorchid4"))
# 
# # Plotting
# geos = subset(geo, time > -max(summ$timeMedian))
# 
# (rttPlot = ggplot(summ, aes(x=-timeMedian, y=median2, colour=clade))+
#     geom_vline(xintercept = geos$time, color="lightgrey") +
#     geom_segment(aes(x=-timeMedian-timeSD/2, xend=-timeMedian+timeSD/2, y=median2, yend=median2), 
#                  color="grey90", linewidth=2, lineend="round")+
#     geom_ribbon(aes(ymin=HPDlowMedian, ymax=HPDuppMedian, fill=colour, color=NULL), alpha=0.1)+
#     geom_line()+
#     # scale_y_log10() + annotation_logticks(sides = 'l')+
#     facet_grid(rate ~ root+div, scales="free")+
#     scale_color_manual(values=as.character(sort(unique(summ$colour)))) +
#     scale_fill_manual(values=as.character(sort(unique(summ$colour)))) +
#     scale_x_continuous(breaks=seq(-max(round(summ$timeMedian, -2)), 0, 100),
#                        minor_breaks=seq(-max(round(summ$timeMedian, -2)), 0, 50)) +
#     labs(y="Rates Through Time", x="Time before present")+
#     theme_classic())
# 
# # Export plot
# pdf(paste0(files$dirPlots, "/", files$RTTplot), width=11.69, height=8.27, paper='special'); plot(rttPlot); dev.off()
# 
# 
# 
# tmp = subset(summ, clade!="Alveolata" & clade!="Holozoa" & clade!="Nucletmycea" & clade!="Stramenopila")
# (rttPlots = ggplot(tmp, aes(x=-timeMedian, y=median2, colour=clade))+
#     geom_vline(xintercept = geos$time, color="lightgrey") +
#     geom_segment(aes(x=-timeMedian-timeSD/2, xend=-timeMedian+timeSD/2, y=median2, yend=median2), 
#                  color="grey90", linewidth=2, lineend="round")+
#     geom_ribbon(aes(ymin=HPDlowMedian, ymax=HPDuppMedian, fill=colour, color=NULL), alpha=0.1)+
#     geom_line()+
#     # scale_y_log10() + annotation_logticks(sides = 'l')+
#     facet_grid(rate ~ root+div, scales="free")+
#     scale_color_manual(values=as.character(sort(unique(tmp$colour)))) +
#     scale_fill_manual(values=as.character(sort(unique(tmp$colour)))) +
#     scale_x_continuous(breaks=seq(-max(round(tmp$timeMedian, -2)), 0, 100),
#                        minor_breaks=seq(-max(round(tmp$timeMedian, -2)), 0, 50)) +
#     labs(y="Rates Through Time", x="Time before present")+
#     theme_classic())
# pdf(paste0(files$dirPlots, "/", sub("\\.pdf", "_subset.pdf", files$RTTplot)), width=11.69, height=8.27, paper='special'); plot(rttPlots); dev.off()
# 
# library(tidyr)
# summl = pivot_wider(summ, id_cols=c(root, clade, div, timePoint, timeMedian), names_from=rate, values_from=median2)
# summl$turnover = summl$extinction / summl$speciation
# 
# tmp = pivot_longer(summl, cols=c(diversification, extinction, speciation, turnover), names_to="rate", values_to="value")
# (rttPlots = ggplot(tmp, aes(x=-timeMedian, y=value, colour=clade))+
#     geom_vline(xintercept = geos$time, color="lightgrey") +
#     geom_line()+
#     # scale_y_log10() + annotation_logticks(sides = 'l')+
#     facet_grid(rate ~ root+div, scales="free")+
#     scale_color_manual(values=as.character(sort(unique(summ$colour)))) +
#     scale_fill_manual(values=as.character(sort(unique(summ$colour)))) +
#     scale_x_continuous(breaks=seq(-max(round(summ$timeMedian, -2)), 0, 100),
#                        minor_breaks=seq(-max(round(summ$timeMedian, -2)), 0, 50)) +
#     labs(y="Rates Through Time", x="Time before present")+
#     theme_classic())
# pdf(paste0(files$dirPlots, "/", sub("\\.pdf", "_turnover.pdf", files$RTTplot)), width=11.69, height=8.27, paper='special'); plot(rttPlots); dev.off()
# 
# 
# #

#----  
#---- Reading files --------------------------------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo", "files", "getBLbin")])

# Reading all shift files
data = data.frame()
for(file in files$shifts){
  cat("\r  Working on file (", grep(file, files$shifts), "/", length(files$shifts), ") ", file, "                    ", sep="")
  tmp = fread(file)
  tmp$file = file
  tmp$calibration = sub("-2\\/.*", "", file)
  tmp$root = file %>% sub("MC0(1|2)-2\\/", "", .) %>% sub("\\/.*", "", .)
  tmp$run = file %>% sub(".*root.\\/", "", .) %>% sub("\\/.*", "", .)
  tmp$clade = file %>% sub(".*diversification_", "", .) %>% sub("_div.*", "", .)
  tmp$diversity = as.numeric(file %>% sub(".*div", "", .) %>% sub("_.*", "", .))
  data = rbind(data, tmp)
}; rm(file, tmp); cat("\rDone", rep(" ", 100), "\n")

table(data$calibration)
table(sub(".*_", "", unique(paste0(data$run, "_", data$clade))))

# Create another variable with minimum and maximum diversity
data$div = NA
for(clad in unique(data$clade)){
  ss = subset(data, clade==clad)
  m = mean(unique(ss$diversity))
  data$div[which(data$clade==clad & data$diversity < m)] = "min"
  data$div[which(data$clade==clad & data$diversity > m)] = "max"
}; rm(clad, ss, m)
# tmp=unique(select(data, c("file", "clade", "diversity", "div")))

# Estimate the shifts on diversification rate
# lam(t) = lam1 * exp(lam2 * t)
# mu(t) = mu1 * exp(mu2 * t)
data$diversification = data$lam1 - data$mu1

# Correcting ultrametric mistakes
if(sum(data$timeRev > 0)){data$timeRev = ifelse(data$timeRev > 0, 0, data$timeRev)}

data$eon = fifelse(data$timeRev < -541, "Proterozoic", "Phanerozoic")
table(data$eon)

data$era = data$timeRev
for(e in rev(unique(geo$era))){
  m = ifelse(e == "Cenozoic", 0, max(subset(geo, era==e)$timeLow))
  data$era[which(data$timeRev <= m)] = e
}
table(data$era)

data$period = data$timeRev
for(e in rev(unique(geo$period))){
  m = ifelse(e == "Quaternary", 0, max(subset(geo, period==e)$timeLow))
  data$period[which(data$timeRev <= m)] = e
}; rm(e, m)
table(data$period)

#
#---- Get whether the shift corresponds to a decay or a growth and correct number of shifts --------

data$change = NA
bl = data.frame(); treenl = c()
for(f in unique(data$file)){
  cat("\r  Working on (", grep(f, unique(data$file)), "/", length(unique(data$file)), ") ", f, "                    ", sep="")
  dirRoot = f %>% sub("\\/bamm/.*", "", .)
  clade = f %>% sub(".*diversification_", "", .) %>% sub("_div.*", "", .)
  ss = subset(data, file==f)
  treen = files$trees %>% grep(dirRoot, ., value=TRUE) %>% grep(clade, ., value=TRUE)
  tree = read.tree(treen) 
  # Get branch length
  calibration = sub("-2\\/.*", "", treen)
  root = treen %>% sub("MC0(1|2)-2\\/", "", .) %>% sub("\\/.*", "", .)
  run = treen %>% sub(".*root.\\/", "", .) %>% sub("\\/clades.*", "", .)
  if(!treen %in% treenl){
    treenl = c(treenl, treen)
    for(g in unique(c(geo$eon, geo$era, geo$period))){
      bl = rbind(bl, data.frame(calibration=calibration, root=root, run=run, clade=clade, time=g, bl=getBLbin(tree, geoBin=g)))
    }
  }
  # Now check whether the shift corresponds to growth or decay
  for(n in 1:length(ss$node)){
    tmpp = paste0("\t", round(n/length(ss$node)*100), "%          ")
    cat("\r  Working on (", grep(f, unique(data$file)), "/", length(unique(data$file)), ") ", f, tmpp, sep="")
    if(ss$time[n] == 0){
      tmp = "root"
    }else{
      noden = ss$node[n]
      dn = subset(ss, node==noden)$diversification
      ancestors = Ancestors(tree, noden)
      if(any(ancestors %in% ss$node)){
        parent = max(ancestors[ancestors %in% ss$node])
        dp = subset(ss, node==parent)$diversification
        tmp = ifelse(dn > dp, "growth", "decay")
      }else{
        tmp = NA
      }
    }
    data$change[which(data$file==f & data$index == n)] = tmp
  }; cat("\r")
}; rm(f, dirRoot, clade, ss, treen, treenl, tree, root, run, g, n, noden, dn, ancestors, parent, dp, tmp, tmpp); cat("\nDone\n")


table(data$change)
if(any(is.na(data$change))){
  cat("\n\n\n########################################\n\n  There is a problem with clades",
      unique(subset(data, is.na(change))$clade), "in", 
      unique(subset(data, is.na(change))$run), "!!\n\n########################################\n\n\n")
}


#
#---- Writing and reading files --------------------------------------------------------------------

# write.table(data, paste0(files$dirPlots, "/", files$shiftsData), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
# write.table(bl, paste0(files$dirPlots, "/", sub("\\.tsv", "_bl.tsv", files$shiftsData)), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

data = fread(paste0(files$dirPlots, "/", files$shiftsData))
bl = fread(paste0(files$dirPlots, "/", sub("\\.tsv", "_bl.tsv", files$shiftsData)))

# Correct diversity categories
data$div2 = NA
for(clad in unique(data$clade)){
  if(clad=="Alveolata"){
    m = 50
  }else if(clad=="Amoebozoa"){
    m = 78
  }else if(clad=="Archaeplastida"){
    m = 56
  }else if(clad=="Cryptista"){
    m = 28
  }else if(clad=="Discoba"){
    m = 68
  }else if(clad=="Haptista"){
    m = 65
  }else if(clad=="Holozoa"){
    m = 70
  }else if(clad=="Metamonada"){
    m = 65
  }else if(clad=="Nucletmycea"){
    m = 75
  }else if(clad=="Rhizaria"){
    m = 50
  }else if(clad=="Stramenopila"){
    m = 50
  }
  data$div2[which(data$clade==clad & data$diversity < m)] = "min"
  data$div2[which(data$clade==clad & data$diversity > m)] = "max"
}; rm(clad, m)
sum(data$div == data$div2)
nrow(data)
data$div = data$div2; data$div2 = NULL

data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
data$change = factor(data$change, levels=c("root", "growth", "decay"))
data$calibration = factor(data$calibration, levels=c("MC01", "MC02"))
data$root = factor(data$root, levels=c("rootD", "rootA"))

#
# #---- Get Branch lengths ---------------------------------------------------------------------------

# tree = read.tree("rootD/fRAcat1-1/clades/clade_Rhizaria.tre")
# getBLbin(tree, geoBin="Ediacaran", verbose=TRUE)
# getBLbin(tree, geoBin="Cambrian", verbose=TRUE)
# getBLbin(tree, geoBin="Quaternary", verbose=TRUE)
# # for(i in c("Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary")){getBLbin(tree, geoBin=i, verbose=TRUE)};rm(i)
# getBLbin(tree, 1000, 541, verbose=TRUE)
# rm(tree)

# bl = data.frame()
# for(f in files$trees){
#   cat("\r  Working on (", grep(f, unique(files$trees)), "/", length(unique(files$trees)), ") ", f, "                    ", sep="")
#   tree = read.tree(f)
#   root = sub("\\/.*", "", f)
#   run = f %>% sub("root.\\/", "", .) %>% sub("\\/clades.*", "", .)
#   clade = f %>% sub(".*clade_", "", .) %>% sub("\\.tre", "", .)
#   for(g in unique(c(geo$eon, geo$era, geo$period))){
#     bl = rbind(bl, data.frame(root=root, run=run, clade=clade, time=g, bl=getBLbin(tree, geoBin=g)))
#   }
#   if(FALSE){
#     seqdates = seq(max(node.depth.edgelength(tree)), 0, length.out=51)
#     for(g in 1:50){
#       timen = ifelse(g<10, paste0("i0", g), paste0("i", g))
#       bl = rbind(bl, data.frame(root=root, run=run, clade=clade, time=timen, bl=getBLbin(tree, max=seqdates[g], min=seqdates[g+1])))
#     }
#   }
# }; rm(f, tree, root, run, clade, g, seqdates, timen); cat("\nDone\n")

data$eonBL = NA; data$eraBL = NA; data$periodBL = NA
for(i in 1:nrow(bl)){
  cat("\r  ", round(i/nrow(bl)*100, 1), "%    ", sep="")
  if(bl$time[i] %in% data$eon){
    data$eonBL[which(data$clade==bl$clade[i] & data$root==bl$root[i] & data$run==bl$run[i] & data$eon==bl$time[i])] = bl$bl[i]
  }else if(bl$time[i] %in% data$era){
    data$eraBL[which(data$clade==bl$clade[i] & data$root==bl$root[i] & data$run==bl$run[i] & data$era==bl$time[i])] = bl$bl[i]
  }else if(bl$time[i] %in% data$period){
    data$periodBL[which(data$clade==bl$clade[i] & data$root==bl$root[i] & data$run==bl$run[i] & data$period==bl$time[i])] = bl$bl[i]
  }
};cat("\r  Done\n"); rm(i)

write.table(data, paste0(files$dirPlots, "/", files$shiftsData), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(bl, paste0(files$dirPlots, "/", files$bl), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# data = fread(paste0(files$dirPlots, "/", files$shiftsData))
# bl = fread(paste0(files$dirPlots, "/", files$bl))

#----  
#---- Checking number of shifts and average branch lengths per period ------------------------------

# Plot the avergae number of shifts per period
datas = data
datas$period = ifelse(datas$period=="Quaternary" | datas$period=="Neogene" | datas$period=="Paleogene", "Cenozoic", datas$period)

datas = subset(datas, change !="root") %>% group_by(calibration, root, clade, run, div, period, change) %>% summarise(shifts=length(unique(index)))
datas = datas %>% group_by(calibration, root, clade, div, period, change) %>% 
  summarise(replicateS=length(unique(run)), sm=mean(shifts, na.rm=TRUE), smsd=sd(shifts, na.rm=TRUE), st=sum(shifts))

unique(datas$period)
datas$period = factor(datas$period, levels=c("Orosirian", "Statherian", 
                                             "Calymnian", "Ectasian", "Stenian", 
                                             "Tonian", "Cryogenian", "Ediacaran", 
                                             "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", 
                                             "Triassic", "Jurassic", "Cretaceous", "Cenozoic"))
datas$clade = factor(datas$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa",
                                           "Metamonada","Discoba",
                                           "Haptista","Cryptista","Archaeplastida",
                                           "Rhizaria","Stramenopila","Alveolata"))

{datas$clades = as.character(datas$clade)
  datas$clades[which(datas$clade=="Amoebozoa")] = "Am"
  datas$clades[which(datas$clade=="Nucletmycea")] = "N"
  datas$clades[which(datas$clade=="Holozoa")] = "Ho"
  datas$clades[which(datas$clade=="Metamonada")] = "M"
  datas$clades[which(datas$clade=="Discoba")] = "D"
  datas$clades[which(datas$clade=="Haptista")] = "Ha"
  datas$clades[which(datas$clade=="Cryptista")] = "C"
  datas$clades[which(datas$clade=="Archaeplastida")] = "Ar"
  datas$clades[which(datas$clade=="Rhizaria")] = "R"
  datas$clades[which(datas$clade=="Stramenopila")] = "S"
  datas$clades[which(datas$clade=="Alveolata")] = "Al"
  datas$clades = factor(datas$clades, levels=c("Am","N","Ho","M","D","Ha","C","Ar","R","S","Al"))}

{datas$colour = as.character(datas$clade)
  datas$colour[which(datas$colour=="Amoebozoa")]="royalblue1"
  datas$colour[which(datas$colour=="Nucletmycea")]="steelblue2"
  datas$colour[which(datas$colour=="Holozoa")]="steelblue4"
  datas$colour[which(datas$colour=="Metamonada")]="forestgreen"
  datas$colour[which(datas$colour=="Discoba")]="orange2"
  datas$colour[which(datas$colour=="Haptista")]="yellow1"
  datas$colour[which(datas$colour=="Cryptista")]="hotpink1"
  datas$colour[which(datas$colour=="Archaeplastida")]="darkseagreen3"
  datas$colour[which(datas$colour=="Rhizaria")]="darkorchid2"
  datas$colour[which(datas$colour=="Stramenopila")]="darkorchid3"
  datas$colour[which(datas$colour=="Alveolata")]="darkorchid4"
  datas$colour = factor(datas$colour,
                        levels=c("royalblue1","steelblue2","steelblue4",
                                 "forestgreen","orange2",
                                 "yellow1","hotpink1","darkseagreen3",
                                 "darkorchid2", "darkorchid3","darkorchid4"))}


data$calibration = factor(data$calibration, levels=c("MC01", "MC02"))
datas$root = factor(datas$root, levels=c("rootD", "rootA"))
datas$change = factor(datas$change, levels=c("growth", "decay"))

datas$st = ifelse(datas$st==1, 1.5, datas$st)

(shiftsNumber = ggplot(datas)+
    # geom_bar(aes(x=clades, y=st, fill=clade), stat="identity", position=position_dodge())+
    geom_bar(data=subset(datas, div=="min"), aes(x=clades, y=st, fill=clade), stat="identity", position=position_dodge())+
    geom_point(data=subset(datas, div=="max"), aes(x=clades, y=st, colour=clade))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    # facet_grid(root+div+change~period, scales="free")+
    facet_grid(root+calibration+change~period, scales="free")+
    scale_fill_manual(values=as.character(sort(unique(datas$colour)))) +
    scale_colour_manual(values=as.character(sort(unique(datas$colour)))) +
    theme_minimal()+
    theme(axis.text.x=element_text(angle=30, hjust=1),
          legend.position = "none"))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_shift_numbers.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsNumber); dev.off()


# Plot the average branch length per period
dataBL = bl
dataBL = subset(dataBL, time %in% geo$period)
dataBL$period = ifelse(dataBL$time=="Quaternary" | dataBL$time=="Neogene" | dataBL$time=="Paleogene", "Cenozoic", dataBL$time)
dataBL = subset(dataBL, !is.na(bl))

dataBL = dataBL %>% group_by(calibration, root, clade, period) %>% summarise(replicateBL=length(unique(run)), blt=sum(bl, na.rm=TRUE), 
                                                                blm=mean(bl, na.rm=TRUE), sd=sd(bl, na.rm=TRUE))
# dataBL = subset(dataBL, replicateBL > 25)

unique(dataBL$period)
dataBL$period = factor(dataBL$period, levels=c("Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Cenozoic"))
(shiftsBL = ggplot(dataBL)+
    geom_point(aes(x=clade, y=blt, alpha=0.8))+
    # geom_point(aes(x=clade, y=blm, alpha=0.8))+
    # geom_errorbar(aes(x=clade, ymin=blm-(sd/2), ymax=blm+(sd/2)),
    #               stat="identity", width=0, position=position_dodge())+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    facet_grid(root+calibration~period, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=c("#FFB000", "#648FFF"))+
    theme(axis.text.x=element_text(angle=30, hjust=1)))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_shift_branchLengths.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsBL); dev.off()

rm(datas, shiftsNumber, dataBL, shiftsBL)

#
#----  
#---- Get corrected number of shifts by time periods -----------------------------------------------

tmp1 = data
tmp1 = subset(tmp1, time != 0)
tmp1$period = ifelse(tmp1$period=="Quaternary" | tmp1$period=="Neogene" | tmp1$period=="Paleogene", "Cenozoic", tmp1$period)
tmp1 = tmp1 %>% group_by(calibration, root, run, clade, div, change, period) %>%
  summarise(shifts=length(unique(index)))

tmp2 = bl
tmp2 = subset(tmp2, time %in% tmp1$period)
tmp2$bl = ifelse(is.na(tmp2$bl), 0, tmp2$bl)
colnames(tmp2) = c("calibration", "root", "run", "clade", "period", "bl")
tmp3 = tmp2
tmp3$div="min"
tmp2$div="max"
tmp2 = rbind(tmp2, tmp3);rm(tmp3)

datac = merge(tmp1, tmp2, by=c("calibration", "root", "run", "clade", "period", "div"), all=TRUE); rm(tmp1, tmp2)
datac = datac %>% group_by(calibration, root, clade, div, period, change) %>% 
  summarise(sm=mean(shifts, na.rm=TRUE), ssd=sd(shifts, na.rm=TRUE), 
            blm=mean(bl, na.rm=TRUE), blsd=sd(bl, na.rm=TRUE), 
            st=sum(shifts, na.rm=TRUE), blt=sum(bl, na.rm=TRUE),
            rep=length(unique(run)))

# datac$scm = datac$sm / datac$blm
datac$sct = datac$st / datac$blt
# datac$scm2 = ifelse(datac$change=="decay", -datac$scm, datac$scm)
datac$sct2 = ifelse(datac$change=="decay", -datac$sct, datac$sct)
datac$sct2t = ifelse(datac$change=="decay", -sqrt(datac$sct), sqrt(datac$sct))

{datac$clades = as.character(datac$clade)
  datac$clades[which(datac$clade=="Amoebozoa")] = "Am"
  datac$clades[which(datac$clade=="Nucletmycea")] = "N"
  datac$clades[which(datac$clade=="Holozoa")] = "Ho"
  datac$clades[which(datac$clade=="Metamonada")] = "M"
  datac$clades[which(datac$clade=="Discoba")] = "D"
  datac$clades[which(datac$clade=="Haptista")] = "Ha"
  datac$clades[which(datac$clade=="Cryptista")] = "C"
  datac$clades[which(datac$clade=="Archaeplastida")] = "Ar"
  datac$clades[which(datac$clade=="Rhizaria")] = "R"
  datac$clades[which(datac$clade=="Stramenopila")] = "S"
  datac$clades[which(datac$clade=="Alveolata")] = "Al"
  datac$clades = factor(datac$clades, levels=c("Am","N","Ho","M","D","Ha","C","Ar","R","S","Al"))}

{datac$colour = as.character(datac$clade)
  datac$colour[which(datac$colour=="Amoebozoa")]="royalblue1"
  datac$colour[which(datac$colour=="Nucletmycea")]="steelblue2"
  datac$colour[which(datac$colour=="Holozoa")]="steelblue4"
  datac$colour[which(datac$colour=="Metamonada")]="forestgreen"
  datac$colour[which(datac$colour=="Discoba")]="orange2"
  datac$colour[which(datac$colour=="Haptista")]="yellow1"
  datac$colour[which(datac$colour=="Cryptista")]="hotpink1"
  datac$colour[which(datac$colour=="Archaeplastida")]="darkseagreen3"
  datac$colour[which(datac$colour=="Rhizaria")]="darkorchid2"
  datac$colour[which(datac$colour=="Stramenopila")]="darkorchid3"
  datac$colour[which(datac$colour=="Alveolata")]="darkorchid4"
  datac$colour = factor(datac$colour,
                        levels=c("royalblue1","steelblue2","steelblue4",
                                 "forestgreen","orange2",
                                 "yellow1","hotpink1","darkseagreen3",
                                 "darkorchid2", "darkorchid3","darkorchid4"))}

# ggplot(datac, aes(y=rep, x=0))+geom_boxplot()+geom_violin(alpha=0.4)+geom_jitter(alpha=0.2)
# datac = subset(datac, replicateS > 2 & replicateBL > 2)

datac$period = factor(datac$period, levels=c("Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Cenozoic"))

(shiftsCperiodS = ggplot(subset(datac, rep>0), aes(x=clades, y=sct2, fill=clade, colour=clade))+
    # geom_bar(stat="identity", position="stack")+
    geom_bar(data=subset(datac, div=="min"), stat="identity", position="stack")+
    geom_point(data=subset(datac, div=="max"))+
    facet_grid(root+calibration~period, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=as.character(sort(unique(datac$colour)))) +
    scale_colour_manual(values=as.character(sort(unique(datac$colour)))) +
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=6),
          legend.position = "none"))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_periods_stacked.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsCperiodS); dev.off()


(shiftsCperiodS = ggplot(subset(datac, rep>0), aes(x=clades, y=sct2t, fill=clade, colour=clade))+
    # geom_bar(stat="identity", position="stack")+
    geom_bar(data=subset(datac, div=="min"), stat="identity", position="stack")+
    geom_point(data=subset(datac, div=="max"))+
    facet_grid(root+calibration~period, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=as.character(sort(unique(datac$colour)))) +
    scale_colour_manual(values=as.character(sort(unique(datac$colour)))) +
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=6),
          legend.position = "none"))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_periods_stacked_sqrt.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsCperiodS); dev.off()

tmp = datac %>% group_by(calibration, period) %>% summarise(sctm=mean(sct, na.rm=TRUE), sctsd=sd(sct, na.rm=TRUE), sctq50=median(sct, na.rm=TRUE))
tmp$period = factor(tmp$period, levels=c("Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacaran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Cenozoic"))
(shiftsCperiodM = ggplot(tmp, aes(x=period, y=sctm))+
    geom_point()+
    geom_errorbar(aes(ymin=sctm-sctsd/2, ymax=sctm+sctsd/2))+
    facet_grid(calibration~., scales="free")+
    theme_bw())
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_periods_stacked_sqrt_mean.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsCperiodM); dev.off()

#
#---- Get corrected number of shifts by eras -------------------------------------------------------

tmp1 = data
tmp1 = subset(tmp1, time != 0)
tmp1 = tmp1 %>% group_by(root, run, clade, div, change, era) %>%
  summarise(shifts=length(unique(index)))

tmp2 = bl
tmp2 = subset(tmp2, time %in% tmp1$era)
tmp2$bl = ifelse(is.na(tmp2$bl), 0, tmp2$bl)
colnames(tmp2) = c("root", "run", "clade", "era", "bl")
tmp3 = tmp2
tmp3$div="min"
tmp2$div="max"
tmp2 = rbind(tmp2, tmp3);rm(tmp3)

datac = merge(tmp1, tmp2, by=c("root", "run", "clade", "era", "div"), all=TRUE); rm(tmp1, tmp2)
datac = datac %>% group_by(root, clade, div, era, change) %>% 
  summarise(sm=mean(shifts, na.rm=TRUE), ssd=sd(shifts, na.rm=TRUE), 
            blm=mean(bl, na.rm=TRUE), blsd=sd(bl, na.rm=TRUE), 
            st=sum(shifts, na.rm=TRUE), blt=sum(bl, na.rm=TRUE),
            rep=length(unique(run)))

# datac$scm = datac$sm / datac$blm
datac$sct = datac$st / datac$blt
# datac$scm2 = ifelse(datac$change=="decay", -datac$scm, datac$scm)
datac$sct2 = ifelse(datac$change=="decay", -datac$sct, datac$sct)
datac$sct2t = ifelse(datac$change=="decay", -sqrt(datac$sct), sqrt(datac$sct))

{datac$clades = as.character(datac$clade)
  datac$clades[which(datac$clade=="Amoebozoa")] = "Am"
  datac$clades[which(datac$clade=="Nucletmycea")] = "N"
  datac$clades[which(datac$clade=="Holozoa")] = "Ho"
  datac$clades[which(datac$clade=="Metamonada")] = "M"
  datac$clades[which(datac$clade=="Discoba")] = "D"
  datac$clades[which(datac$clade=="Haptista")] = "Ha"
  datac$clades[which(datac$clade=="Cryptista")] = "C"
  datac$clades[which(datac$clade=="Archaeplastida")] = "Ar"
  datac$clades[which(datac$clade=="Rhizaria")] = "R"
  datac$clades[which(datac$clade=="Stramenopila")] = "S"
  datac$clades[which(datac$clade=="Alveolata")] = "Al"
  datac$clades = factor(datac$clades, levels=c("Am","N","Ho","M","D","Ha","C","Ar","R","S","Al"))}

{datac$colour = as.character(datac$clade)
  datac$colour[which(datac$colour=="Amoebozoa")]="royalblue1"
  datac$colour[which(datac$colour=="Nucletmycea")]="steelblue2"
  datac$colour[which(datac$colour=="Holozoa")]="steelblue4"
  datac$colour[which(datac$colour=="Metamonada")]="forestgreen"
  datac$colour[which(datac$colour=="Discoba")]="orange2"
  datac$colour[which(datac$colour=="Haptista")]="yellow1"
  datac$colour[which(datac$colour=="Cryptista")]="hotpink1"
  datac$colour[which(datac$colour=="Archaeplastida")]="darkseagreen3"
  datac$colour[which(datac$colour=="Rhizaria")]="darkorchid2"
  datac$colour[which(datac$colour=="Stramenopila")]="darkorchid3"
  datac$colour[which(datac$colour=="Alveolata")]="darkorchid4"
  datac$colour = factor(datac$colour,
                        levels=c("royalblue1","steelblue2","steelblue4",
                                 "forestgreen","orange2",
                                 "yellow1","hotpink1","darkseagreen3",
                                 "darkorchid2", "darkorchid3","darkorchid4"))}

# ggplot(datac, aes(y=rep, x=0))+geom_boxplot()+geom_violin(alpha=0.4)+geom_jitter(alpha=0.2)
# datac = subset(datac, replicateS > 2 & replicateBL > 2)

datac$era = factor(datac$era, levels=c("Mesoproterozoic", "Neoproterozoic", "Paleozoic", "Mesozoic", "Cenozoic"))

(shiftsCeraS = ggplot(subset(datac, rep>0), aes(x=clades, y=sct2, fill=clade))+
    geom_bar(stat="identity", position="stack")+
    facet_grid(root+div~era, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=as.character(sort(unique(datac$colour)))) +
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=6)))
# pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_eras_stacked.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsCperiodS); dev.off()


(shiftsCeraS = ggplot(subset(datac, rep>0), aes(x=clades, y=sct2t, fill=clade))+
    geom_bar(stat="identity", position="stack")+
    facet_grid(root+div~era, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=as.character(sort(unique(datac$colour)))) +
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=6)))
# pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_eras_stacked_sqrt.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsCperiodS); dev.off()

tmp = datac %>% group_by(era) %>% summarise(sctm=mean(sct, na.rm=TRUE), sctsd=sd(sct, na.rm=TRUE), sctq50=median(sct, na.rm=TRUE))
tmp$era = factor(tmp$era, levels=c("Mesoproterozoic", "Neoproterozoic", "Paleozoic", "Mesozoic", "Cenozoic"))
ggplot(tmp, aes(x=era, y=sctm))+geom_point()+geom_errorbar(aes(ymin=sctm-sctsd/2, ymax=sctm+sctsd/2))+theme_minimal()

#
#----
#---- Checking number of shifts and average branch lengths per eons --------------------------------


tmp1 = data
tmp1 = subset(tmp1, change !="root") %>% group_by(root, clade, run, div, eon, change) %>% summarise(shifts=length(unique(index)))
tmp1 = tmp1 %>% group_by(root, clade, div, eon, change) %>% summarise(replicateS=length(unique(run)), shiftsM=mean(shifts, na.rm=TRUE), shiftsSD=sd(shifts, na.rm=TRUE), shiftsT=sum(shifts))

tmp2 = bl
tmp2 = subset(tmp2, time %in% geo$eon)
tmp2$root = gsub("\\/.*", "", tmp2$tree)
tmp2 = subset(tmp2, !is.na(bl))
tmp2 = tmp2 %>% group_by(root, clade, time) %>% summarise(replicateBL=length(unique(tree)), blm=mean(bl, na.rm=TRUE), sd=sd(bl, na.rm=TRUE))
colnames(tmp2) = c("root", "clade", "eon", "replicateBL", "blm", "sd")
tmp = merge(tmp1, tmp2, by=c("root", "clade", "eon"))

tmp$shiftsC = tmp$shiftsM/tmp$blm
tmp$shiftsCt = tmp$shiftsT/tmp$blm
tmp$shiftsCt2 = ifelse(tmp$change == "decay", -tmp$shiftsCt, tmp$shiftsCt)

tmp$eon = factor(tmp$eon, levels=c("Proterozoic", "Phanerozoic"))

(shiftsNumberCperiodS = ggplot(tmp, aes(x=clade, y=shiftsCt2, fill=change))+
    geom_bar(stat="identity", position="stack")+
    facet_grid(root+div~eon, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=c("#648FFF", "#FFB000"))+
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=8)))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_eons_stacked.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsNumberCperiodS); dev.off()

#
#---- Checking number of shifts and average branch lengths per eras --------------------------------

tmp1 = data
tmp1 = subset(tmp1, change !="root") %>% group_by(root, clade, run, div, era, change) %>% summarise(shifts=length(unique(index)))
tmp1 = tmp1 %>% group_by(root, clade, div, era, change) %>% summarise(replicateS=length(unique(run)), shiftsM=mean(shifts, na.rm=TRUE), shiftsSD=sd(shifts, na.rm=TRUE), shiftsT=sum(shifts))

tmp2 = bl
tmp2 = subset(tmp2, time %in% geo$era)
tmp2$root = gsub("\\/.*", "", tmp2$tree)
tmp2 = subset(tmp2, !is.na(bl))
tmp2 = tmp2 %>% group_by(root, clade, time) %>% summarise(replicateBL=length(unique(tree)), blm=mean(bl, na.rm=TRUE), sd=sd(bl, na.rm=TRUE))
colnames(tmp2) = c("root", "clade", "era", "replicateBL", "blm", "sd")

tmp = merge(tmp1, tmp2, by=c("root", "clade", "era")); rm(tmp1, tmp2)

tmp$shiftsC = tmp$shiftsM/tmp$blm
tmp$shiftsCt = tmp$shiftsT/tmp$blm
tmp$shiftsCt2 = ifelse(tmp$change == "decay", -tmp$shiftsCt, tmp$shiftsCt)

tmp$era = factor(tmp$era, levels=c("Mesoproterozoic", "Neoproterozoic", "Paleozoic", "Mesozoic", "Cenozoic"))

(shiftsNumberCperiodS = ggplot(tmp, aes(x=clade, y=shiftsCt2, fill=change))+
    geom_bar(stat="identity", position="stack")+
    facet_grid(root+div~era, scales="free")+
    theme_minimal()+
    scale_fill_manual(values=c("#648FFF", "#FFB000"))+
    # theme(axis.text.x=element_text(angle=30, hjust=1, size=8)))
    theme(axis.text.x=element_text(size=8)))
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_corrected_eras_stacked.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsNumberCperiodS); dev.off()

#
#----
#----
#---- Plot shifts as bulk --------------------------------------------------------------------------

if(!"data" %in% ls()){data = fread(paste0(files$dirPlots, "/", gsub("\\.pdf$", "_rawData.tsv", files$shiftsPlot)))}

geos = subset(geo, time>min(data$timeRev))

# Add factors
data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
data$shift = ifelse(data$time == 0, "root", "shift")
data$shift = factor(data$shift, levels=c("root", "shift"))
data$change = factor(data$change, levels=c("root", "growth", "decay"))

# plot
(shifts = ggplot(data, aes(x=timeRev, y=clade, fill=shift))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_violin(scale="width", alpha=0.6)+
    # geom_point(aes(colour=change), alpha=0.1)+
    facet_wrap(root~div, ncol=2)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title="Shifts in diversification rates",
         y="Clades", x="Time before present (Ma)")+
    scale_fill_manual(values=c("grey20", "turquoise4")) +
    # scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic())

pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allCombined_noScenario.pdf"), width=11.69, height=8.27, paper='special'); plot(shifts); dev.off()
#

# plot
data = subset(data, !is.na(change))
(shiftsPulledChange = ggplot(data, aes(x=timeRev, y=change, fill=change))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    # geom_violin(scale="area", alpha=0.6)+
    # geom_violin(scale="count", alpha=0.6)+
    geom_violin(scale="width", alpha=0.6)+
    # geom_point(aes(colour=change), alpha=0.1)+
    facet_wrap(root~div, ncol=2)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title="Shifts in diversification rates",
         y="Clades", x="Time before present (Ma)")+
    scale_fill_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    # scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic())

pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allCombined_change.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsPulledChange); dev.off()

# plot
(shiftsPulled = ggplot(data, aes(x=timeRev, y=shift, fill=shift))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_violin(scale="width", alpha=0.6)+
    # geom_point(aes(colour=change), alpha=0.1)+
    facet_wrap(root~div, ncol=2)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title="Shifts in diversification rates",
         y="Clades", x="Time before present (Ma)")+
    scale_fill_manual(values=c("grey20", "turquoise4")) +
    # scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic())

pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allCombined.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsPulled); dev.off()
#

#---- And plot growth or decays --------------------------------------------------------------------

if(!"data" %in% ls()){data = fread(paste0(files$dirPlots, "/", gsub("\\.pdf$", "_rawData.tsv", files$shiftsPlot)))}

geos = subset(geo, time>min(data$timeRev))

# Add factors
data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
data$root = factor(data$root, levels=c("rootD", "rootA"))
data$change = factor(data$change, levels=c("root", "growth", "decay"))

# plot
datas = subset(data, !is.na(change) & change != "root")
(shifts = ggplot(datas, aes(x=timeRev, y=clade, fill=change))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_violin(scale="width", alpha=0.6)+
    # geom_point(aes(colour=change), alpha=0.1)+
    # facet_wrap(div~root, ncol=2)+
    facet_grid(div+root~.)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title="Shifts in diversification rates",
         y="Clades", x="Time before present (Ma)")+
    scale_fill_manual(values=c("#648FFF", "#FFB000")) +
    # scale_fill_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    # scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic()+
    theme(legend.position = "none"))

pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allShifts.pdf"), height=11.69*1.5, width=8.27, paper='special'); plot(shifts); dev.off()

# Now only shifts happening in the mesozoic or older, only for Rhizaria and Holozoa and with at least 90% of the replicates
# timeLimit = -250
timeLimit = -541
# timeLimit = 0

# datas = subset(data, timeRev <= timeLimit & (clade == "Rhizaria" | clade == "Holozoa"))
datas = subset(data, timeRev <= timeLimit & !is.na(change))
geos = subset(geo, time <= timeLimit & time>min(datas$timeRev))

(shiftsSubset = ggplot(datas, aes(x=timeRev, y=clade, fill=change))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_violin(scale="width", alpha=0.4)+
    # geom_point(aes(size=lam1, colour=change), alpha=0.2)+
    facet_wrap(~div, ncol=1)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title="Shifts in diversification rates on the Proterozoic",
         y="Clades", x="Time before present (Ma)")+
    scale_fill_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic())
# geom_violin; scale: if "area" (default), all violins have the same area (before trimming the tails). If "count", areas are scaled proportionally to the number of observations. If "width", all violins have the same maximum width.
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allShifts_proterozoic.pdf"), width=11.69, height=8.27, paper='special'); plot(shiftsSubset); dev.off()

# And finally, one file at a time
data$run = data$file %>% sub("\\/[^\\/]+$", "", .) %>% sub("\\/bamm/.*", "", .)
pdf(paste0(files$dirPlots, "/", files$shiftsPlot, "_allRuns.pdf"), width=11.69, height=8.27, paper='special')
for(r in unique(data$run)){
  cat("\r  Plotting ", "(", grep(r, unique(data$run)), "/", length(unique(data$run)), ") ", r, "                    ", sep="")
  tmp = subset(data, run == r)
  geos = subset(geo, time>-max(tmp$time))
  # tmp = subset(tmp, timeRev < -251)
  # geos = subset(geos, time < -251)
  shiftsR = ggplot(tmp, aes(x=timeRev, y=clade, size=diversification, colour=change))+
    geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_point(alpha=0.4)+
    facet_wrap(~div, ncol=1)+
    scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                       minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
    labs(title=paste0("Shifts in diversification rates: run ", r),
         y="Clades", x="Time before present (Ma)")+
    scale_colour_manual(values=c("grey20", "#648FFF", "#FFB000")) +
    theme_classic()
  plot(shiftsR)
}; rm(r, tmp); cat("\rDone                                        \n"); dev.off()

#

#----  
#---- Others ---------------------------------------------------------------------------------------

library(geiger)
library(data.table)
library(dplyr)

setwd("/home/miguel/Documents/Uppsala/1_mecoEvo/data/euk/stepDating/MC01/rootD/")

diversity = fread("fRAcat1-1/clads/fractions.tsv")

tree = read.tree("fRAcat1-1/clades/clade_Archaeplastida.tre")
div = max(subset(diversity, V1=="clade_Archaeplastida.tre")$V2)

richness = data.frame(taxon=tree$tip.label, 
                      n.taxa=rep(round((1-div)+1), length(tree$tip.label)))
head(richness)

res = medusa(tree, richness, verbose=TRUE,
             model="bd",
             init = c(r=0.0009, epsilon=0.3),
             ncores=8)

plot(res, cex=0.5, label.offset=1, edge.width=2) 


#----  