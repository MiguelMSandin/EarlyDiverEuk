#! /bin/Rscript

# loading packages _________________________________________________________________________________

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(TESS))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(ggplot2))

# Set parser _______________________________________________________________________________________
parser = OptionParser()

parser = add_option(parser, c("-t", "--tree"), dest="tree", type="character",
                    help="The time-tree in newick format.")

parser = add_option(parser, c("-f", "--fraction"), dest="fraction", type="double", default=1,
                    help="The sampling fraction [0-1], default='1'.")

parser = add_option(parser, c("-m", "--mcmc"), dest="mcmc", type="integer", default=200000,
                    help="The maximum number of iteration of the MCMC, default='200000'.")

parser = add_option(parser, c("-d", "--dir"), dest="dir", type="character", default="NONE",
                    help="The output directory. By default, will remove the extension of the input file.")

parser = add_option(parser, c("-p", "--prefix"), dest="prefix", type="character", default="NONE",
                    help="The prefix file name for the plot and the data. By default, will replace the extension of the input file by '.pdf' and '.tsv.")

parser = add_option(parser, c("-v", "--verbose"), dest="verbose", action="store_true", default=TRUE,
                    help="If selected, will not print information to the console.")

args = parse_args(parser)

# Setting file names _______________________________________________________________________________
if(args$dir== "NONE"){
  dirOut = sub("\\.[^\\.]+$", "", args$tree)
}else{
  dirOut = args$dir
}
if(args$prefix== "NONE"){
  outPlot = sub("\\.[^\\.]+$", ".pdf", args$tree)
  outData = sub("\\.[^\\.]+$", ".tsv", args$tree)
}else{
  outPlot = paste0(args$prefix, ".pdf")
  outData = paste0(args$prefix, ".tsv")
}

geo = data.frame(time= c(2500, 2300, 2050, 1800, 1600, 1400, 1200, 1000, 720, 635, 541, 485.4, 443.8, 419, 358.9, 298.9, 251.9, 201.4, 145, 66, 23.03, 2.58),
                 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
                 era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

era = data.frame(time= c(2500, 1600, 1000, 541, 251.9, 66),
                 era=c("Paleoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Paleozoic", "Mesozoic", "Cenozoic"))

# Reading tree _____________________________________________________________________________________
if(args$verbose){cat("  Reading file:", args$tree, "\n")}
tree = read.tree(args$tree)

if(!is.ultrametric(tree)){
  if(args$verbose){cat("    Ultrametricizing tree\n")}
  tree = force.ultrametric(tree, message=FALSE)
}

# Running CoMET ____________________________________________________________________________________
seed=round(runif(1, 0, 10^9))
if(args$verbose){cat("  Running CoMET\n    Sampling fraction:", args$fraction, "\n    MCMC chain length:", args$mcmc,"\n    Seed:", seed,"\n")}

if(!dir.exists(dirOut)){
  if(args$verbose){cat("    Creating output directory:", dirOut,"\n")}
  dir.create(dirOut, recursive=TRUE)
}

set.seed(seed)
tess.analysis(tree,
              empiricalHyperPriors = TRUE,
              MAX_ITERATIONS = args$mcmc,
              dir = dirOut)

# Saving data ______________________________________________________________________________________
if(args$verbose){cat("  Exporting data to:", outData, "\n")}

output = tess.process.output(dirOut)

data = data.frame(time=head(output$intervals, -1),
                  pp=colMeans(output$`mass extinction times`),
                  bf=output$`mass extinction Bayes factors`)

write.table(data, outData, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Plotting _________________________________________________________________________________________
if(args$verbose){cat("  Exporting plot to:", outPlot, "\n")}

geos = subset(geo, time < max(data$time))
geoss = subset(era, time < max(data$time))

plotME = ggplot(data, aes(x=time, y=pp))+
  geom_hline(yintercept = output$massExtinctionCriticalPosteriorProbabilities, color="grey90", linetype = "dashed") +
  geom_vline(xintercept = geos$time, color="grey80") +
  geom_vline(xintercept = geoss$time, color="grey50") +
  geom_col() + 
  scale_x_reverse(breaks=seq(round(max(data$time), -2), 0, -100),
                     minor_breaks=seq(round(max(data$time), -2), 0, -50)) + 
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(title=paste0(args$tree),
       subtitle=paste0("Sampling fraction:", args$fraction),
       x="Million years ago",
       y="Posterior probability")

pdf(outPlot, width=11.69, height=8.27, paper='special')
plot(plotME)
dev.off()

# Plotting _________________________________________________________________________________________
if(args$verbose){cat("Done\n")}
