#----
#---- loading packages and data.frames -------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, 
                          -1600 , -1400, -1200, 
                          -1000, -720, -635, 
                          -541, -485.4, -443.8, -419, -358.9, -298.9, 
                          -251.9, -201.4, -145, 
                          -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", 
                           "Calymnian","Ectasian","Stenian",
                           "Tonian","Cryogenian","Ediacran",
                           "Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian",
                           "Triassic","Jurassic","Cretaceous",
                           "Paleogene","Neogene","Quaternary"),
                  era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic",
                        "Mesoproterozoic","Mesoproterozoic","Mesoproterozoic",
                        "Neoproterozoic","Neoproterozoic","Neoproterozoic",
                        "Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
                        "Mesozoic","Mesozoic","Mesozoic",
                        "Cenozoic","Cenozoic","Cenozoic"))
geo$mid <- apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean) 

#----
#---- set working directory and file names ---------------------------------------------------------
setwd("~/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01r/")

files = list(dir="rootD/",
             extension="LTT-data.tsv",
             factor=-1,
             out="rootDiscoba_cleaned_LTT_221128") # Without extension for the tsv table and the pdf

rm(list=ls()[!ls() %in% c("geo", "files")])

#----
#---- Read tables ----------------------------------------------------------------------------------

tmp = grep(paste0("*", files$extension), dir(files$dir, recursive=TRUE), value=TRUE)
data = data.frame()
i = 0
for(file in tmp){
    cat("\r  Reading file ", (i = i + 1), "/", length(tmp), sep="")
    f <- fread(paste0(files$dir, file))
    f$file = gsub("/.*", "", file)
    data = rbind(data, f)
}; rm(tmp, i, file, f); cat("\n")

#----
#----
#---- Get stats ------------------------------------------------------------------------------------

# Get slopes per 100 years of every supergroup and replicate

slopes = data.frame()
for(group in unique(data$tree)){
    i = 0
    tmp = data[data$tree == group,]
    if(max(tmp$N >= 300)){
        for(replicate in unique(data$file)){
            cat("\r  Analysing ", group, " (", (i = i + 1), "/", length(unique(data$file)), ")          ", sep="")
            tmp1 = tmp[tmp$file == replicate,]
            M = round(min(tmp1$time))+100
            j = 0
            for(m in c(seq(M, 0, 100), 0)){
                j = j + 100
                tmp2 = subset(tmp1, time <= m)
                slopes <- rbind(slopes, data.frame(group=group, file=replicate, 
                                                   age=m, meassured = j,
                                                   N=max(tmp2$N), logN=max(tmp2$logN),
                                                   slope=lm(tmp2$logN ~ tmp2$time)$coeff[2]))
            }
        }
    }else{
        cat("  ", group, " not analysed: less than 300 lineages in tree.", sep="")
    }
    cat("\n")
}; rm(i, group, replicate, tmp, tmp1, M, m, tmp2, j); cat("Done\n")


sort(unique(slopes$group))
slopes$group <- factor(slopes$group, 
                    levels=c("Main_tree", 
                             "Amoebozoa", "Nucletmycea", "Holozoa", 
                             "Metamonada", "Discoba", 
                             "Haptista", "Cryptista", "Archaeplastida", 
                             "Rhizaria", "Stramenopila", "Alveolata"))

slopes$colour <- as.character(slopes$group)
{slopes$colour[which(slopes$colour=="Main_tree")]="black"
    slopes$colour[which(slopes$colour=="Amoebozoa")]="royalblue1"
    slopes$colour[which(slopes$colour=="Apusomonadida")]="lightskyblue1"
    # slopes$colour[which(slopes$colour=="Breviatea")]="lightskyblue2"
    slopes$colour[which(slopes$colour=="Nucletmycea")]="steelblue2"
    slopes$colour[which(slopes$colour=="Holozoa")]="steelblue4"
    # slopes$colour[which(slopes$colour=="Ancyromonadida")]="grey20"
    # slopes$colour[which(slopes$colour=="CRuMs")]="aquamarine1"
    slopes$colour[which(slopes$colour=="Metamonada")]="forestgreen"
    slopes$colour[which(slopes$colour=="Discoba")]="orange2"
    slopes$colour[which(slopes$colour=="Haptista")]="yellow1"
    slopes$colour[which(slopes$colour=="Cryptista")]="hotpink1"
    # slopes$colour[which(slopes$colour=="Hemimastigophora")]="grey80"
    slopes$colour[which(slopes$colour=="Archaeplastida")]="darkseagreen3"
    # slopes$colour[which(slopes$colour=="Telonemia")]="plum4"
    slopes$colour[which(slopes$colour=="Rhizaria")]="darkorchid2"
    slopes$colour[which(slopes$colour=="Stramenopila")]="darkorchid3"
    slopes$colour[which(slopes$colour=="Alveolata")]="darkorchid4"}

slopes$colour <- factor(slopes$colour, 
                      levels=c("black",
                               "royalblue1", "lightskyblue1", "steelblue2", "steelblue4",
                               "forestgreen", "orange2",
                               "yellow1", "hotpink1", "darkseagreen3",
                               "darkorchid2", "darkorchid3", "darkorchid4"))

slopes$meassured <- as.character(slopes$meassured)

ggplot(slopes, aes(x=age, y=slope, colour=group))+
    # geom_point()+
    geom_jitter(alpha=0.4)+
    geom_smooth()+
    theme_bw()+
    scale_color_manual(values=as.character(sort(unique(slopes$colour))))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

tmp = subset(slopes, meassured == 100 | 
                 meassured == 200 | 
                 meassured == 300 | 
                 meassured == 400)

(slopes = ggplot(tmp, aes(x=group, y=slope, colour=meassured))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))

pdf("LTT_slopes_timeIntervals.pdf", width=11.69, height=8.27, paper='special')
plot(slopes)
dev.off()


#----
#----
#---- Summarise, subset and export table -----------------------------------------------------------

# Summarise table
data$time = data$time * files$factor
summ <- data %>% group_by(tree, N, logN) %>% summarise(# mean=mean(time),
                                                       # sd=sd(time),
                                                       # min=min(time),
                                                       q05=quantile(time, 0.05),
                                                       # q25=quantile(time, 0.25),
                                                       median=quantile(time, 0.50),
                                                       # q75=quantile(time, 0.75),
                                                       q95=quantile(time, 0.95),
                                                       # max=max(time),
                                                       replicates=length(unique(file)))

# only use lineages and time points with at least 50% of the replicates
tmp <- max(summ$replicates)
summ <- subset(summ, replicates >= tmp*0.5)

# Only use groups with at least 300 lineages
summ = summ[summ$tree %in% names(table(summ$tree)[table(summ$tree) > 300]),]

# Exporting the dataset ____________________________________________________________________________
write.table(summ, paste0(files$out, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
#----
#---- Reading table if already available -----------------------------------------------------------

summ = fread(paste0(files$out, ".tsv"))

#----
#----
#---- beautifying the dataset for plotting ---------------------------------------------------------

sort(unique(summ$tree))
summ$tree <- factor(summ$tree, 
                    levels=c("Main_tree", 
                             "Amoebozoa", "Nucletmycea", "Holozoa", 
                             "Metamonada", "Discoba", 
                             "Haptista", "Cryptista", "Archaeplastida", 
                             "Rhizaria", "Stramenopila", "Alveolata"))

summ$colour <- as.character(summ$tree)
{summ$colour[which(summ$colour=="Main_tree")]="black"
    summ$colour[which(summ$colour=="Amoebozoa")]="royalblue1"
    summ$colour[which(summ$colour=="Apusomonadida")]="lightskyblue1"
    # summ$colour[which(summ$colour=="Breviatea")]="lightskyblue2"
    summ$colour[which(summ$colour=="Nucletmycea")]="steelblue2"
    summ$colour[which(summ$colour=="Holozoa")]="steelblue4"
    # summ$colour[which(summ$colour=="Ancyromonadida")]="grey20"
    # summ$colour[which(summ$colour=="CRuMs")]="aquamarine1"
    summ$colour[which(summ$colour=="Metamonada")]="forestgreen"
    summ$colour[which(summ$colour=="Discoba")]="orange2"
    summ$colour[which(summ$colour=="Haptista")]="yellow1"
    summ$colour[which(summ$colour=="Cryptista")]="hotpink1"
    # summ$colour[which(summ$colour=="Hemimastigophora")]="grey80"
    summ$colour[which(summ$colour=="Archaeplastida")]="darkseagreen3"
    # summ$colour[which(summ$colour=="Telonemia")]="plum4"
    summ$colour[which(summ$colour=="Rhizaria")]="darkorchid2"
    summ$colour[which(summ$colour=="Stramenopila")]="darkorchid3"
    summ$colour[which(summ$colour=="Alveolata")]="darkorchid4"}

summ$colour <- factor(summ$colour, 
                      levels=c("black",
                               "royalblue1", "lightskyblue1", "steelblue2", "steelblue4",
                               "forestgreen", "orange2",
                               "yellow1", "hotpink1", "darkseagreen3",
                               "darkorchid2", "darkorchid3", "darkorchid4"))

#----
#---- plotting -------------------------------------------------------------------------------------

# First subset the geological dates table
geos <- subset(geo, time > -max(summ$q95, na.rm=TRUE))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(lttplot <- ggplot(summ)+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_line(aes(x=-median, y=logN, color=tree)) +
        geom_ribbon(aes(xmin=-q95, xmax=-q05, y=logN, fill=tree), alpha=0.2, colour = NA) +
        # geom_line(aes(x=-q05, y=logN, color=tree), alpha=0.2) +
        # geom_line(aes(x=-q95, y=logN, color=tree), alpha=0.2) +
        geom_text(data=geos, aes(x=mid, y=max(summ$logN)+1, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(summ$logN)+0.5, label=period), size=2)+
        scale_x_continuous(breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 100), 
                           minor_breaks=seq((round(min(-summ$q95, na.rm=TRUE), -2)), 0, 50)) +
        scale_y_continuous(breaks=1:(round(max(summ$logN), 0))) +
        scale_color_manual(values=as.character(sort(unique(summ$colour))))+
        scale_fill_manual(values=as.character(sort(unique(summ$colour))))+
        theme_classic()+
        labs(y="Ln (N)", x="Time (Ma)"))

pdf(paste0(files$out, "_summary.pdf"), width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off()

#----
#----
