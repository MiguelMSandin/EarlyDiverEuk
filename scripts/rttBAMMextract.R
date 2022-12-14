#! /bin/Rscript

# loading packages _________________________________________________________________________________

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(BAMMtools))
suppressPackageStartupMessages(library(HDInterval))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# Set parser _______________________________________________________________________________________
parser <- OptionParser()

parser <- add_option(parser, c("-t", "--tree"), dest="tree", type="character",
                     help="Tree file name in nexus format with annotated nodes.")

parser <- add_option(parser, c("-e", "--eventData"), dest="edata", type="character",
                     help="The event data file.")

parser <- add_option(parser, c("-o", "--out"), dest="out", type="character", default="NONE",
                     help="The output tab delimited table. By default, will add '_RTT.tsv' to the input event data file (after removing '_event_data.txt').")

parser <- add_option(parser, c("-b", "--burnin"), dest="burnin", type="integer", default=25,
                     help="The burnin (in percentage; 0-100), default='25'.")

parser <- add_option(parser, c("-d", "--dates"), dest="dates", action="store_false", default=FALSE,
                     help="If selected, will extract the starting dates from the root of the tree and not from the BAMM analyses.")

parser <- add_option(parser, c("-p", "--plot"), dest="plot", action="store_true", default=TRUE,
                     help="If selected, will not print a tree with the coloured branch specific rates to 'outputname.pdf'.")

parser <- add_option(parser, c("-v", "--verbose"), dest="verbose", action="store_true", default=TRUE,
                     help="If selected, will not print information to the console.")

args = parse_args(parser)

# Setting file names _______________________________________________________________________________
if(args$out== "NONE"){
    output = sub("_event_data.txt", "_RTT.tsv", args$edata)
}else{
    output = args$out
}
if(args$plot){
    outplot = sub("\\.[^\\.]+$", "_RTTtree.pdf", output)
}

# Open the analyzed tree ___________________________________________________________________________
if(args$verbose){cat("  Reading tree:", args$tree, "\n")}
tree <- read.tree(args$tree)
# plot(tree, cex=.6)
# axisPhylo()

# Reading event data _______________________________________________________________________________
if(args$verbose){cat("  Reading event data:", args$edata,  "\n    Burnin:", args$burnin, "\n")}
edata <- getEventData(args$tree, 
                      eventdata = args$edata, 
                      burnin=args$burnin/100)

if(args$plot){
    if(args$verbose){cat("  Plotting branch specific rates to:", outplot, "\n")}
    pdf(outplot, width=11.69, height=8.27, paper='special')
    bamm.tree = plot.bammdata(edata, lwd=2, labels=TRUE, cex=0.2, logcolor=TRUE)
    addBAMMshifts(edata, cex=2)
    addBAMMlegend(bamm.tree)
    title(main = args$edata)
    dev.off()
}

# Extracting RTT matrix ____________________________________________________________________________
if(args$dates){
    if(args$verbose){cat("  Extracting RTT matrix \n")}
    edatam <- getRateThroughTimeMatrix(edata)
}else{
    # or set times manually
    start.time = -ltt.plot.coords(tree)[1]
    if(args$verbose){cat("  Extracting RTT matrix from", start.time, "to", 0,"units before present\n")}
    edatam <- getRateThroughTimeMatrix(edata, start.time=start.time, end.time=0)
}



# Create a table for speciation data
tmp <- as.data.frame(t(edatam$lambda))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_spe <- cbind(time=rev(edatam$times), tmp)

# Create a table for extinction data
tmp <- as.data.frame(t(edatam$mu))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_ext <- cbind(time=rev(edatam$times), tmp)

# Create a table for diversification data
tmp <- as.data.frame(t(edatam$lambda)) - as.data.frame(t(edatam$mu))
colnames(tmp) <- paste0("replicate", seq(1:ncol(tmp)))
edatam_div <- cbind(time=rev(edatam$times), tmp)

if(args$verbose){cat("    Summarising data\n")}
# first summarise speciation data
hpd = apply(edatam_spe[,-1], 1, hdi)
edatamOut = data.frame(rate="speciation",
                       time=edatam_spe$time,
                       median=apply(edatam_spe[,-1], 1, median),
                       HPDlow=hpd[1,],
                       HPDupp=hpd[2,],
                       mean=apply(edatam_spe[,-1], 1, mean),
                       sd=apply(edatam_spe[,-1], 1, sd))

# now summarise extinction data
hpd = apply(edatam_ext[,-1], 1, hdi)
edatamOut = rbind(edatamOut, data.frame(rate="extinction",
                                        time=edatam_ext$time,
                                        median=apply(edatam_ext[,-1], 1, median),
                                        HPDlow=hpd[1,],
                                        HPDupp=hpd[2,],
                                        mean=apply(edatam_ext[,-1], 1, mean),
                                        sd=apply(edatam_ext[,-1], 1, sd)))

# and finally summarise diversification data
hpd = apply(edatam_div[,-1], 1, hdi)
edatamOut = rbind(edatamOut, data.frame(rate="diversification",
                                        time=edatam_div$time,
                                        median=apply(edatam_div[,-1], 1, median),
                                        HPDlow=hpd[1,],
                                        HPDupp=hpd[2,],
                                        mean=apply(edatam_div[,-1], 1, mean),
                                        sd=apply(edatam_div[,-1], 1, sd)))

# Extracting RTT matrix ____________________________________________________________________________
if(args$verbose){cat("  Exporting summary file\n")}
write.table(edatamOut, output, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

if(args$verbose){cat("Done\n")}
