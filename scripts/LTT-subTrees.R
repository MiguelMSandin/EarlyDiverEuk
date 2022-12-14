
#---- loading packages and data.frames

library(ape)
suppressMessages(library(treeio))
library(BAMMtools)

library(HDInterval)

library(data.table)

library(ggplot2)

library(optparse)

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

# setwd("/home/miguel/Desktop/Uppsala/1_ecoEvo/data/euk/stepDating/MC01/rootD/RAgamma1-1/")

parser <- OptionParser()

parser <- add_option(parser, c("-t", "--tree"), dest="tree", type="character",
                help="Tree file name in nexus format.")

parser <- add_option(parser, c("-o", "--out"), dest="out", type="character", default="NONE",
                help="The output tab delimited table. By default, will add '_LTT-data.tsv' to the input tree (after removing the extension).")

parser <- add_option(parser, c("-f", "--factor"), dest="factor", type="integer", default=-1,
                help="The factor to multiply the node heights, default='-1'.")

# parser <- add_option(parser, c("-v", "--verbose"), action="store_true", default=TRUE,
#               help="If selected, will not print information to the console")

opt = parse_args(parser)


cat("  Reading tree\n")

if(opt$out== "NONE"){
    output = sub("\\.[^\\.]+$", "_LTT-data.tsv", opt$tree)
}else{
    output = opt$out
}

treei <- read.beast(opt$tree)
tree <- as.phylo(treei)

# Transform files
treed <- as_tibble(treei)
names(treed) <- gsub("!", "", names(treed))

# Select all annotated nodes _______________________________________________________________________
annotations <- sort(as.vector(treed$name[!is.na(treed$name)]))


annotations[!annotations %in% c("Amoebozoa", "Apusomonadida", "Breviatea", "Nucletmycea", "Holozoa", "Ancyromonadida", "CRuMs", "Metamonada", "Discoba", "Haptista", "Cryptista", "Hemimastigophora", "Archaeplastida", "Telonemia", "Rhizaria", "Stramenopila", "Alveolata")]

# Loop through the main tree and the different annotated nodes to get the LTT data _________________
cat("  Extracting LTT data\n")
lttData <- data.frame()
for(annot in c("Main_tree", annotations)){
    cat("\r    Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")                    ", sep="")
    if(annot=="Main_tree"){
        tmp <- as.data.frame(subset(treed, is.na(label)))
    }else{
        node <- as.numeric(subset(treed, name==annot)$node)
        tmp <- tree_subset(treei, node=node, levels_back=0)
        tmp <- as_tibble(tmp)
        tmp <- as.data.frame(subset(tmp, is.na(label)))
    }
    tmp <- tmp[order(as.numeric(tmp$height_median), decreasing = TRUE),]
    tmp$N <- seq(1, nrow(tmp))
    tmp <- data.frame(tree=annot,
                      time=as.numeric(tmp$height_median),
                      N=tmp$N,
                      hpd=as.character(tmp$height_0.95_HPD))
    lttData <- rbind(lttData, tmp)
}
rm(annot, tmp, node)
cat("\r    Done                                   \n")

# Transform variables ______________________________________________________________________________
cat("  Transforming variables\n")
lttData$logN <- log(lttData$N, base=exp(1))
lttData$time <- lttData$time * opt$factor
lttData <- subset(lttData, !is.na(time))
lttData$hpd_min <- ifelse(lttData$tree=="Main_tree", as.numeric(lttData$hpd %>% gsub("c\\(", "", .) %>% gsub(",.*", "", .)) * opt$factor, NA)
lttData$hpd_max <- ifelse(lttData$tree=="Main_tree", as.numeric(lttData$hpd %>% gsub("\\)", "", .) %>% gsub(".*, ", "", .)) * opt$factor, NA)

# Load dataset if already available________________________________________________________________
# lttData <- fread(gsub("\\.tre$", "_LTT-data.tsv", files$tree))

# Get summary per annotated groups _________________________________________________________________
lineages <- data.frame()
for(annot in unique(lttData$tree)){
    ss <- subset(lttData, tree==annot)
    lineages <- rbind(lineages, data.frame(group=annot,
                                           lineages=max(ss$N),
                                           origin=min(ss$time),
                                           slope=lm(time ~ N, ss)$coefficients[2],
                                           slopeLog=lm(time ~ logN, ss)$coefficients[2]))
}; rm(annot, ss)

# Check if all annotations are in the factorization ________________________________________________
factorization <- c("Main_tree", "Amoebozoa", "Apusomonadida", "Breviatea", "Nucletmycea", "Holozoa", "Ancyromonadida", "CRuMs", "Metamonada", "Discoba", "Haptista", "Cryptista", "Hemimastigophora", "Archaeplastida", "Telonemia", "Rhizaria", "Stramenopila", "Alveolata")
if(!all(annotations %in% factorization)){annotations[!annotations %in% factorization]}

if(!all(factorization[-1] %in% annotations)){factorization[!factorization %in% annotations]}

# Beautifying the dataset for the plot _____________________________________________________________
sort(unique(lttData$tree))
lttData$tree <- factor(lttData$tree, 
                       levels=c("Main_tree", 
                                "Amoebozoa", "Apusomonadida", "Breviatea", "Nucletmycea", "Holozoa", 
                                "Ancyromonadida", "CRuMs", "Metamonada", "Discoba", 
                                "Haptista", "Cryptista", "Hemimastigophora", "Archaeplastida", 
                                "Telonemia", "Rhizaria", "Stramenopila", "Alveolata"))
c("black",
  "royalblue1", "lightskyblue1", "lightskyblue2", "steelblue2", "steelblue4",
  "grey20", "aquamarine1", "forestgreen", "orange2",
  "yellow1", "hotpink1", "grey80", "darkseagreen3",
  "plum4", "darkorchid2", "darkorchid3", "darkorchid4")

lttData$colour <- as.character(lttData$tree)
{lttData$colour[which(lttData$colour=="Main_tree")]="black"
    lttData$colour[which(lttData$colour=="Amoebozoa")]="royalblue1"
    lttData$colour[which(lttData$colour=="Apusomonadida")]="lightskyblue1"
    lttData$colour[which(lttData$colour=="Breviatea")]="lightskyblue2"
    lttData$colour[which(lttData$colour=="Nucletmycea")]="steelblue2"
    lttData$colour[which(lttData$colour=="Holozoa")]="steelblue4"
    lttData$colour[which(lttData$colour=="Ancyromonadida")]="grey20"
    lttData$colour[which(lttData$colour=="CRuMs")]="aquamarine1"
    lttData$colour[which(lttData$colour=="Metamonada")]="forestgreen"
    lttData$colour[which(lttData$colour=="Discoba")]="orange2"
    lttData$colour[which(lttData$colour=="Haptista")]="yellow1"
    lttData$colour[which(lttData$colour=="Cryptista")]="hotpink1"
    lttData$colour[which(lttData$colour=="Hemimastigophora")]="grey80"
    lttData$colour[which(lttData$colour=="Archaeplastida")]="darkseagreen3"
    lttData$colour[which(lttData$colour=="Telonemia")]="plum4"
    lttData$colour[which(lttData$colour=="Rhizaria")]="darkorchid2"
    lttData$colour[which(lttData$colour=="Stramenopila")]="darkorchid3"
    lttData$colour[which(lttData$colour=="Alveolata")]="darkorchid4"}

lttData$colour <- factor(lttData$colour, 
                         levels=c("black",
                                  "royalblue1", "lightskyblue1", "lightskyblue2", "steelblue2", "steelblue4",
                                  "grey20", "aquamarine1", "forestgreen", "orange2",
                                  "yellow1", "hotpink1", "grey80", "darkseagreen3",
                                  "plum4", "darkorchid2", "darkorchid3", "darkorchid4"))


# Exporting the dataset ____________________________________________________________________________
cat("  Exporting table\n")
write.table(lttData, output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Plot the raw dataset _____________________________________________________________________________

data <- lttData

geos <- subset(geo, time > min(c(data$hpd_max, data$time), na.rm=TRUE))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

# And plotting _____________________________________________________________________________________
cat("  Plotting\n")
lttplot <- ggplot(data, aes(x=time, y=logN))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_segment(aes(x=hpd_min, xend=hpd_max, y=logN, yend=logN), color="grey90", size=2, lineend="round")+
        geom_line(aes(color=tree)) +
        # geom_point(aes(color=tree)) +
        geom_text(data=geos, aes(x=mid, y=max(data$logN)+1, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(data$logN)+0.5, label=period), size=2)+
        scale_x_continuous(breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 100), 
                           minor_breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 50)) +
        scale_y_continuous(breaks=1:(round(max(data$logN), 0))) +
        scale_color_manual(values=as.character(sort(unique(data$colour))))+
        theme_classic()+
        labs(y="Ln (N)", x="Time (Ma)")

pdf(gsub("\\.tre$", "_LTT-subtress.pdf", opt$tree), width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off()

# Selecting the dataset to be plotted ______________________________________________________________

data <- lttData[with(lttData, tree %in% subset(lineages, lineages>300)$group),]
# data <- subset(lttData, tree=="Main_tree" | tree=="Rhizaria" | tree=="Holozoa")

geos <- subset(geo, time > min(c(data$hpd_max, data$time), na.rm=TRUE))
geos$mid <- apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

# And plotting _____________________________________________________________________________________
lttplot <- ggplot(data, aes(x=time, y=logN))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_segment(aes(x=hpd_min, xend=hpd_max, y=logN, yend=logN), color="grey90", size=2, lineend="round")+
        geom_line(aes(color=tree)) +
        # geom_point(aes(color=tree)) +
        geom_text(data=geos, aes(x=mid, y=max(data$logN)+1, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(data$logN)+0.5, label=period), size=2)+
        scale_x_continuous(breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 100), 
                           minor_breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 50)) +
        scale_y_continuous(breaks=1:(round(max(data$logN), 0))) +
        scale_color_manual(values=as.character(sort(unique(data$colour))))+
        theme_classic()+
        labs(y="Ln (N)", x="Time (Ma)")

pdf(gsub("\\.tre$", "_LTT-subtressClean.pdf", opt$tree), width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off() 

cat("Done\n")
