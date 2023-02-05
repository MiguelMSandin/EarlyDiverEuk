#---- 
#---- loading packages ----
# 
# library(ape)
# library(phangorn)
# library(treeio)
# 
# library(BAMMtools)
# library(HDInterval)
# library(coda)
# 
# library(data.table)
# library(tidyr)
# library(dplyr)
# 
# library(ggplot2)
# library(ggtree)
# 
# # library(devtools)
# # install_github("hmorlon/PANDA", dependencies = TRUE)
# library(RPANDA)
# 
# # library(TreePar)
# 
# library(vegan)
# # library(drc)
# # library(nlme)
# # library(aomisc)

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid <- apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----  
setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/repository/9_diversificationCLaDS/clads/")
rm(list=ls()[!ls() %in% c("geo")])
# .rs.restartR()
#---- 
#----
#---- Read Rdata from several processed ClaDS tables (julia) ---------------------------------------

# Checking output from two different ClaDS Julia output ____________________________________________
# .rs.restartR()
rm(list=ls()[!ls() %in% c("geo")])

files = list(all=grep("clads.*RTT\\.tsv", dir(recursive=TRUE), value=TRUE),
			 diversity="max")
files$clades = unique(files$all %>% sub("clads\\/", "", .) %>% sub("\\/.*", "", .))

files$files = c()
for(clade in files$clades){
	f = grep(clade, files$all, value=TRUE)
	if(length(f) == 1){cat("Warning!", clade, "has only one file:", f, "\n")}
	div = f %>% sub(".*_div", "", .) %>% sub("_.*", "", .)
	div = as.numeric(div)
	if(files$diversity == "min"){
		d = min(div)
	} else if(files$diversity == "max"){
		d = max(div)
	}
	files$files = c(files$files, grep(d, f, value=TRUE))
};rm(clade, f, div, d)

files

# files = list(files=c("clads/Amoebozoa/clade_Amoebozoa_div89_clads_RTT.tsv",
# 					  "clads/Nucletmycea/clade_Nucletmycea_div84_clads_RTT.tsv",
# 					  "clads/Holozoa/clade_Holozoa_div83_clads_RTT.tsv",
# 					  "clads/Discoba/clade_Discoba_div84_clads_RTT.tsv",
# 					  "clads/Metamonada/clade_Metamonada_div81_clads_RTT.tsv",
# 					  "clads/Haptista/clade_Haptista_div80_clads_RTT.tsv",
# 					  "clads/Cryptista/clade_Cryptista_div44_clads_RTT.tsv",
# 					  "clads/Archaeplastida/clade_Archaeplastida_div70_clads_RTT.tsv",
# 					  "clads/Rhizaria/clade_Rhizaria_div71_clads_RTT.tsv",
# 					  "clads/Stramenopila/clade_Stramenopila_div70_clads_RTT.tsv",
# 					  "clads/Alveolata/clade_Alveolata_div68_clads.Rdata"),
# 			 diversity="min")

rtt = data.frame()
dtt = data.frame(); DTT=TRUE
for(file in files$files){
    cat("Loading file '", file ,"' (", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    
    tmp = fread(file)
    tmp$revTime = -rev(tmp$time)
    tmp$group = file %>% sub("clads\\/", "", .) %>% sub("\\/.*", "", .)
    rtt = rbind(rtt, tmp)
    
    if(DTT){
        tmp = fread(gsub("RTT", "DTT", file))
        tmp$group = file %>% sub("clads\\/", "", .) %>% sub("\\/.*", "", .)
        tmp$revTime = -rev(tmp$time)
        dtt = rbind(dtt, tmp)
    }
}; rm(tmp, DTT)

geos <- subset(geo, time>-max(rtt$time))

{rtt$colour = rtt$group
    rtt$colour[which(rtt$colour=="Main_tree")]="black"
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

rtt$group = factor(rtt$group, levels=c(#"All",
    "Amoebozoa", 
    "Nucletmycea",
    "Holozoa",
    "Metamonada",
    "Discoba",
    "Haptista",
    "Cryptista",
    "Archaeplastida",
    "Rhizaria",
    "Stramenopila",
    "Alveolata"
    ))

rtt$colour <- factor(rtt$colour,
                     levels=c(#"black",
                         "royalblue1",
                         "steelblue2",
                         "steelblue4",
                         "forestgreen",
                         "orange2",
                         "yellow1",
                         "hotpink1",
                         "darkseagreen3",
                         "darkorchid2", 
                         "darkorchid3",
                         "darkorchid4"
                         ))

(rttplot <- ggplot(rtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_ribbon(aes(x=revTime, ymax=HPD_95, ymin=HPD_05, fill=group), alpha=0.2, colour = NA) +
        geom_line(aes(x=revTime, y=rate, colour=group), size=1)+
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.01, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(rtt$HPD_95)+0.005, label=period), size=2)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(rtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(rtt$colour))))+
        labs(y="Diversification rate (Ln)", x="Time (Ma)"))

pdf(paste0("clads/clads_RTT_div", files$diversity, ".pdf"), width=11.69, height=8.27, paper='special'); plot(rttplot); dev.off()


{dtt$colour = dtt$group
    dtt$colour[which(dtt$colour=="Main_tree")]="black"
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

dtt$group = factor(dtt$group, levels=c(#"All",
    "Amoebozoa",
    "Nucletmycea",
    "Holozoa",
    "Metamonada",
    "Discoba",
    "Haptista",
    "Cryptista",
    "Archaeplastida",
    "Rhizaria",
    "Stramenopila",
    "Alveolata"
))

dtt$colour <- factor(dtt$colour,
                     levels=c(#"black",
                         "royalblue1",
                         "steelblue2",
                         "steelblue4",
                         "forestgreen",
                         "orange2",
                         "yellow1",
                         "hotpink1",
                         "darkseagreen3",
                         "darkorchid2",
                         "darkorchid3",
                         "darkorchid4"
                     ))

(dttplot <- ggplot(dtt)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_ribbon(aes(x=time, ymax=max, ymin=min, fill=group), alpha=0.2, colour = NA) +
        geom_line(aes(x=time, y=lineages, colour=group), size=1)+
        scale_y_log10() + annotation_logticks(sides = 'l')+
        theme_classic()+
        scale_color_manual(values=as.character(sort(unique(dtt$colour))))+
        scale_fill_manual(values=as.character(sort(unique(dtt$colour))))+
        labs(y="Estimated diversity through time (Ln)", x="Time (Ma)"))

pdf(paste0("clads/clads_DTT_div", files$diversity, ".pdf"), width=11.69, height=8.27, paper='special'); plot(dttplot); dev.off()

#
#----
#---- Plot the slope with the increment change on the RTT slope ------------------------------------

files = list(file="clads/Holozoa/clade_Holozoa_div83_clads_RTT.tsv",
             outPlot="clads/Holozoa/changes_RTT_Holozoa.pdf")

# Read all files, and estimate the change in diversification rate
data = fread(files$file)
tmp1 = c(); tmp2 = c()
for(i in 2:nrow(data)){
    tmp2 = c(tmp2, (data$rate[i]-data$rate[(i-1)]))
}; rm(i)
data$change = c(0, tmp2); rm(tmp2)
data$timeRev = -rev(data$time)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))

# Add a colour whether the change is positive or negative
data$colour = ifelse(data$changeSQRT < 0, "negative", "positive")
data$colour = factor(data$colour, levels=c("negative", "positive"))
data$colourSlope = c(data$colour[-1], data$colour[1])
data$colourSlope = factor(data$colourSlope, levels=c("negative", "positive"))

# And now plot
geos <- subset(geo, time>-max(data$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(barPlot = ggplot(data)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                  size=2, lineend = "round")+
        geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                 stat="identity", position=position_dodge(), alpha=0.8)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
        theme_classic()+
        scale_color_manual(values=c("orangered4", "springgreen4")) +
        scale_fill_manual(values=c("orangered4", "springgreen4")) +
        theme(legend.position="none")+
        labs(title=files$file,
             y="Change in diversification rate (square-rooted)", x="Time (Ma)"))

pdf(files$outPlot, width=11.69, height=8.27, paper='special')
plot(barPlot)
dev.off()
#

#---- Plot the slope with the increment change on the RTT slope of different files -----------------

files = list(filesDir="clads/",
             outPlot="clads/changes_RTT.pdf")
files$files = grep("RTT\\.tsv", dir(files$filesDir, recursive=TRUE), value=TRUE)

# Read all files, estimate the change in diversification rate and plot
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
for(file in files$files){
    cat("  Working on '", file, "' (", grep(file, files$files), "/", length(files$files), ")\n", sep="")
    data = fread(paste0(files$filesDir, "/", file))
    tmp1 = c()
    for(i in 2:nrow(data)){
        tmp1 = c(tmp1, (data$rate[i]-data$rate[(i-1)]))
    }
    data$change = c(0, tmp1)
    data$timeRev = -rev(data$time)
    data$file = file
    
    # Square root the change of rates and keep whether it is positive or negative
    data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))
    
    # Add a colour whether the change is positive or negative
    data$colour = ifelse(data$changeSQRT > 0, "positive", ifelse(data$changeSQRT == 0, "zero", "negative"))
    data$colour = factor(data$colour, levels=c("positive", "zero", "negative"))
    data$colourSlope = c(data$colour[-1], data$colour[1])
    data$colourSlope = factor(data$colourSlope, levels=c("positive", "zero", "negative"))
    
    # And now plot
    geos <- subset(geo, time>-max(data$time))
    geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
    
    (barPlot = ggplot(data)+
            geom_vline(xintercept = geos$time, color="lightgrey") +
            geom_line(aes(x=timeRev, y=sqrt(rate), color=colourSlope, group=1), 
                      size=2, lineend = "round")+
            geom_bar(aes(x=timeRev, y=changeSQRT, fill=colour, group=1), 
                     stat="identity", position=position_dodge(), alpha=0.8)+
            # scale_y_log10() + annotation_logticks(sides = 'l')+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
            # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
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
#---- Trying to find diversification-dependent patterns --------------------------------------------

files = list(filesDir="clads/",
             outPlot="clads/changes_RTT_allVsAll.pdf")

files$files = grep("RTT\\.tsv", dir(files$filesDir, recursive=TRUE), value=TRUE)

# Read all files, and estimate the change in diversification rate
data = data.frame(); tmp1 = c(); tmp2 = c()
for(file in files$files){
    tmp = fread(paste0(files$filesDir, "/", file))
    tmp1 = c(); tmp2 = c()
    for(i in 2:nrow(tmp)){
        tmp2 = c(tmp2, (tmp$rate[i]-tmp$rate[(i-1)]))
    }
    tmp$file = file
    tmp$change = c(0, tmp2)
    tmp$timeRev = -rev(tmp$time)
    data = rbind(data, tmp)
}; rm(file, tmp, tmp1, tmp2, i)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))

# And now plot all pair-wise comparisons
pdf(files$outPlot, width=11.69, height=8.27, paper='special')
done = c()
for(file1 in files$files){
    for( file2 in files$files){
        if(file1 != file2 & ((!paste0(file1, "-", file2) %in% done) | (!paste0(file2, "-", file1) %in% done))){
            file1="Rhizaria/clade_Rhizaria_div71_clads_RTT.tsv"
            # file2="Holozoa/clade_Holozoa_div83_clads_RTT.tsv"
            file2="Stramenopila/clade_Stramenopila_div70_clads_RTT.tsv"
            cat("  Plotting", file1, "against", file2, "\n")
            done = c(done, paste0(file1, "-", file2))
            datas = subset(data, file==file1 | file==file2)
            
            geos <- subset(geo, time>-max(datas$time))
            geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)
            
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="pearson")
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="kendall")
            # cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
            #     method="spearman")
            
            (barPlot = ggplot(datas, aes(x=timeMidRev, y=changeSQRT, fill=file, group=file))+
                    geom_vline(xintercept = geos$time, color="lightgrey") +
                    geom_line(aes(x=timeRev, y=sqrt(rate), color=file))+
                    geom_bar(stat="identity", position=position_dodge(), width=10, alpha=0.8)+
                    # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
                    # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
                    theme_classic()+
                    theme(legend.position="bottom")+
                    labs(y="Change in diversification rate (square-rooted)", x="Time (Ma)"))
            plot(barPlot)
        }
    }
}; rm(done, file1, file2, tmp); dev.off()

#

#---- Trying to find diversification-dependent patterns in targetted trees -------------------------

files = list(files=c("clads/Holozoa/clade_Holozoa_div83_clads_RTT.tsv",
                     "clads/Rhizaria/clade_Rhizaria_div71_clads_RTT.tsv"),
             out="clads/changes_RTT_HoloRhiz.pdf")


# Read all files, and estimate the change in diversification rate
data = data.frame(); tmp1 = c(); tmp2 = c()
for(file in files$files){
    tmp = fread(file)
    tmp1 = c(); tmp2 = c()
    for(i in 2:nrow(tmp)){
        tmp2 = c(tmp2, (tmp$rate[i]-tmp$rate[(i-1)]))
    }
    tmp$file = file
    tmp$change = c(0, tmp2)
    tmp$timeRev = -rev(tmp$time)
    data = rbind(data, tmp)
}; rm(file, tmp, tmp1, tmp2, i)

# Square root the change of rates and keep whether it is positive or negative
data$changeSQRT = ifelse(data$change < 0, -sqrt(-data$change), sqrt(data$change))


# Colour the clades

data$clade = data$file %>% sub(".*clade_", "", .) %>% sub("_div.*", "", .)

{data$colour = data$clade
    data$colour[which(data$colour=="Main_tree")]="black"
    data$colour[which(data$colour=="Amoebozoa")]="royalblue1"
    data$colour[which(data$colour=="Nucletmycea")]="steelblue2"
    data$colour[which(data$colour=="Holozoa")]="steelblue4"
    data$colour[which(data$colour=="Metamonada")]="forestgreen"
    data$colour[which(data$colour=="Discoba")]="orange2"
    data$colour[which(data$colour=="Haptista")]="yellow1"
    data$colour[which(data$colour=="Cryptista")]="hotpink1"
    data$colour[which(data$colour=="Archaeplastida")]="darkseagreen3"
    data$colour[which(data$colour=="Rhizaria")]="darkorchid2"
    data$colour[which(data$colour=="Stramenopila")]="darkorchid3"
    data$colour[which(data$colour=="Alveolata")]="darkorchid4"}

# And now plot all pair-wise comparisons
geos <- subset(geo, time>-max(data$time))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

# cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
#     method="pearson")
# cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
#     method="kendall")
# cor(subset(datas, file==file1)$changeSQRT, subset(datas, file==file2)$changeSQRT, 
#     method="spearman")

(barPlot = ggplot(data, aes(x=timeRev, y=changeSQRT, fill=colour, group=clade))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=timeRev, y=sqrt(rate), color=file))+
        # geom_bar(stat="identity", position=position_dodge(), width=10, alpha=0.8)+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.01, label=era), size=2.5)+
        # geom_text(data=geos, aes(x=mid, y=max(datas$changeSQRT)+0.005, label=period), size=2)+
        theme_classic()+
        theme(legend.position="bottom")+
        scale_color_manual(values=as.character(sort(unique(data$colour))))+
        scale_fill_manual(values=as.character(sort(unique(data$colour))))+
        labs(y="Change in diversification rate (square-rooted)", x="Time (Ma)"))
plot(barPlot)


#

#----
# 
