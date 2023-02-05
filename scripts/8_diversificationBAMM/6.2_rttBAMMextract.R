#! /bin/Rscript

# loading packages _________________________________________________________________________________

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(BAMMtools))
suppressPackageStartupMessages(library(ggplot2))
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

parser <- add_option(parser, c("-s", "--shifts"), dest="shifts", type="integer", default=1,
                     help="The expected number of shifts, default='1'.")

parser <- add_option(parser, c("-v", "--verbose"), dest="verbose", action="store_true", default=TRUE,
                     help="If selected, will not print information to the console.")

args = parse_args(parser)

# Setting file names _______________________________________________________________________________
if(args$out== "NONE"){
    output = args$edata %>% sub("\\.gz$", "", .) %>% sub("_event_data.txt", "_RTT.tsv", .)
    outputBest = args$edata %>% sub("\\.gz$", "", .) %>% sub("_event_data.txt", "_RTTbest.tsv", .)
    outputShifts =  args$edata %>% sub("\\.gz$", "", .) %>% sub("_event_data.txt", "_shifts.tsv", .)
}else{
    output = args$out
    outputBest = sub("\\.[^\\.]+$", "_RTTbest.tsv", output)
    outputShifts = sub("\\.[^\\.]+$", "_shifts.tsv", output)
}

# Open the analyzed tree ___________________________________________________________________________
if(args$verbose){cat("  Reading tree:", args$tree, "\n")}
tree <- read.tree(args$tree)

# Reading event data _______________________________________________________________________________
if(args$verbose){cat("  Reading event data:", args$edata,  "\n    Burnin:", args$burnin, "\n")}
edata <- getEventData(args$tree, 
                      eventdata = args$edata, 
                      burnin=args$burnin/100)

# Get best shift configuration _____________________________________________________________________
if(args$verbose){cat("  Getting best shift configuration for", args$shifts, "expected shifts\n")}
best <- getBestShiftConfiguration(edata,
                                  expectedNumberOfShifts=args$shifts,
                                  threshold=5)

shifts = as.data.frame(best$eventData)
shifts$timeRev = shifts$time - max(node.depth.edgelength(tree))
if(args$verbose){cat("  Exporting best shift configuration table to:", outputShifts, "\n")}
write.table(shifts, outputShifts, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# Extracting RTT matrix ____________________________________________________________________________
if(args$verbose){cat("  Extracting RTT matrix \n")}
edatam <- getRateThroughTimeMatrix(edata)

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

# And now export the slopes for the best shift configuration _______________________________________
if(args$verbose){cat("  Extracting RTT matrix for best shift configuration \n")}
edatamb <- getRateThroughTimeMatrix(best)

edatambOut = data.frame(time=rev(edatam$times),
                        speciation=c(edatamb$lambda),
                        extinction=c(edatamb$mu),
                        diversification=c(edatamb$lambda)-c(edatamb$mu))

if(args$verbose){cat("  Exporting summary file\n")}
write.table(edatambOut, outputBest, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

if(args$verbose){cat("Done\n")}
