#---- 
#---- loading packages ----

library(data.table)
library(dplyr)
library(ggplot2)

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid <- apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#---- 
# 'setwd' to the 'bamm/diver' directory of a dated eukaryotic tree, where all the event_data files from different subclades are stored in their respective folders
setwd("")
#---- 
#---- Set file names -------------------------------------------------------------------------------

files = list(files=grep("RTT\\.tsv$", dir(recursive=TRUE), value=TRUE),
			 shifts=grep("shifts\\.tsv$", dir(recursive=TRUE), value=TRUE), 
			 dirTrees="../../clades/",
			 dirPlots="../../plots/",
			 RTTplot="bamm_RTT.pdf",
			 incrementPlot="bamm_RTT_increment.pdf",
			 shiftsPlot="bamm_shifts.pdf")

if(!dir.exists(files$dirPlots)){dir.create(files$dirPlots)}

#---- 
#---- plot RTT BAMM output -------------------------------------------------------------------------

# Reading all RTT files
data = data.frame()
for(file in files$files){
	tmp = fread(file)
	tmp$file = file
	tmp$clade = file %>% sub(".*diversification_", "", .) %>% sub("_.*", "", .)
	tmp$diversity = as.numeric(file %>% sub(".*_div", "", .) %>% sub("_.*", "", .))
	data = rbind(data, tmp)
};rm(file, tmp)

# Create another variable for minimum and maximum diversity
tmp = c()
for(cl in unique(data$clade)){
	ss = subset(data, clade==cl)
	tmp = c(tmp, ifelse(ss$diversity==min(ss$diversity), "min", "max"))
}; rm(cl, ss)
data$div = tmp; rm(tmp)

# Adding colours
{data$colour = data$clade
	data$colour[which(data$colour=="Amoebozoa")] =      "royalblue1"
	data$colour[which(data$colour=="Nucletmycea")] =    "steelblue2"
	data$colour[which(data$colour=="Holozoa")] =        "steelblue4"
	data$colour[which(data$colour=="Metamonada")] =     "forestgreen"
	data$colour[which(data$colour=="Discoba")] =        "orange2"
	data$colour[which(data$colour=="Haptista")] =       "yellow1"
	data$colour[which(data$colour=="Cryptista")] =      "hotpink1"
	data$colour[which(data$colour=="Archaeplastida")] = "darkseagreen3"
	data$colour[which(data$colour=="Rhizaria")] =       "darkorchid2"
	data$colour[which(data$colour=="Stramenopila")] =   "darkorchid3"
	data$colour[which(data$colour=="Alveolata")] =      "darkorchid4"}

# Adding factors for plotting
data$rate = factor(data$rate, levels=c("extinction","speciation","diversification"))
data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))
data$colour = factor(data$colour, levels=c("royalblue1","steelblue2","steelblue4","forestgreen","orange2","yellow1","hotpink1","darkseagreen3","darkorchid2","darkorchid3","darkorchid4"))

# Plotting
(rttPlot = ggplot(data, aes(x=-time, y=median, colour=clade))+
        geom_line()+
        geom_ribbon(aes(ymin=HPDlow, ymax=HPDupp, fill=colour, color=NULL), alpha=0.1)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
		facet_grid(rate ~ div, scales="free")+
		scale_color_manual(values=as.character(sort(unique(data$colour)))) +
		scale_fill_manual(values=as.character(sort(unique(data$colour)))) +
        scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
                           minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
        labs(y="Rates Through Time", x="Time before present")+
        theme_bw())

# Export plot
pdf(paste0(files$dirPlots, "/", files$RTTplot), width=11.69, height=8.27, paper='special'); plot(rttPlot); dev.off()

#

#----  
#---- Plot the slope with the increment change on the RTT slope of different files -----------------

# Read all files, estimate the change in diversification rate and plot
pdf(paste0(files$dirPlots, "/", files$incrementPlot), width=11.69, height=8.27, paper='special')
for(file in files$files){
    cat("  Working on '", file, "' (", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    data = fread(file)
    data = subset(data, rate=="diversification")
    tmp1 = c()
    for(i in 2:nrow(data)){
        tmp1 = c(tmp1, (data$median[i]-data$median[(i-1)]))
    }
    data$change = c(0, tmp1)
    data$timeRev = -data$time
    data$file = file
    
    # Square root the change of rates and keep whether it is positive or negative
    data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))
    
    # Add a colour whether the change is positive or negative
    data$colour = ifelse(data$changeSQRT > 0, "growth", ifelse(data$changeSQRT == 0, "zero", "decay"))
    data$colour = factor(data$colour, levels=c("growth", "zero", "decay"))
    data$colourSlope = c(data$colour[-1], data$colour[1])
    data$colourSlope = factor(data$colourSlope, levels=c("growth", "zero", "decay"))
    
    # And now plot
    geos <- subset(geo, time>-max(data$time))
    geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
    
    (barPlot = ggplot(data)+
            geom_vline(xintercept = geos$time, color="lightgrey") +
            geom_line(aes(x=timeRev, y=sqrt(median), color=colourSlope, group=1), 
            		  linewidth=2, lineend = "round")+
            geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                     stat="identity", position=position_dodge(), alpha=0.8)+
            theme_classic()+
            scale_color_manual(values=c("springgreen4", "black", "orangered4")) +
            scale_fill_manual(values=c("springgreen4", "black", "orangered4")) +
            theme(legend.position="none")+
            labs(title=file,
                 y="Increment in diversification rate (square-rooted)", x="Time (Ma)"))
    plot(barPlot)
}; rm(file, data, tmp1, i); dev.off()


#

#----
#---- Analyze when diversification shifts tend to happen -------------------------------------------

files$trees = grep("clade_.*\\.tre", dir(files$dirTrees), value=TRUE)
files

# Read datasets
data = data.frame()
for(file in files$shifts){
	clade = file %>% sub(".*diversification_", "", .) %>% sub("_div.*", "", .)
	div = file %>% sub(".*div", "", .) %>% sub("_.*", "", .)
	tmp = fread(file)
	tmp$clade = clade
	tmp$diversity = div
	tmp$file = file
	data = rbind(data, tmp)
}; rm(file, clade, div, tmp)

# Create another variable with minimum and maximum diversity
tmp = c()
for(clad in unique(data$clade)){
	ss = subset(data, clade==clad)
	tmp = c(tmp, ifelse(ss$diversity==min(ss$diversity), "min", "max"))
}; rm(clad, ss)
data$div = tmp; rm(tmp)

# Estimate the shifts on diversification rate
# lam(t) = lam1 * exp(lam2 * t)
# mu(t) = mu1 * exp(mu2 * t)
data$diversification = data$lam1 - data$mu2

# Get whether the shift corresponds to a decay or a growth
change = c()
for(f in unique(data$file)){
	cat("\r  Working on (", grep(f, unique(data$file)), "/", length(unique(data$file)), ") ", f, "                    ", sep="")
	ss = subset(data, file==f)
	tree = grep(unique(ss$clade), files$trees, value=TRUE)
	tree = read.tree(paste0(files$dirTrees, "/", tree))
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
};rm(f, ss, tree, changeF, n, noden, dn, ancestors, parent, dp, tmp); cat("\n")

data$change = change; rm(change)
data$change = factor(data$change, levels=c("root", "growth", "decay"))

# And plot
geos <- subset(geo, time>-max(data$time))

# Add factors
data$clade = factor(data$clade, levels=c("Amoebozoa","Nucletmycea","Holozoa","Metamonada","Discoba","Haptista","Cryptista","Archaeplastida","Rhizaria","Stramenopila","Alveolata"))

(shifts = ggplot(data, aes(x=timeRev, y=clade, size=diversification, colour=change))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_point(alpha=0.4)+
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
(shifts = ggplot(subset(data, timeRev < -541), aes(x=timeRev, y=clade, size=diversification, colour=change))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_point(alpha=0.4)+
		facet_wrap(~div, ncol=1)+
		scale_x_continuous(breaks=seq(-max(round(data$time, -2)), 0, 100),
						   minor_breaks=seq(-max(round(data$time, -2)), 0, 50)) +
		labs(title="Shifts in diversification rates on the Proterozoic",
			 y="Clades", x="Time before present (Ma)")+
		scale_colour_manual(values=c("grey20", "springgreen4", "orangered4")) +
		theme_classic())

#

#----  