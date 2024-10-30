#! /bin/Rscript

# loading packages _________________________________________________________________________________
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(geiger))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

# Set parser _______________________________________________________________________________________
parser = OptionParser(description="Given a tree, estimates evolutionary diversification shifts acounting for incomplete lineage sampling and unknown diversity, using the MEDUSA function from the 'geiger' package.")

parser = add_option(parser, c("-t", "--tree"), dest="tree", type="character",
					help="A newick ultrametric tree.")

parser = add_option(parser, c("-d", "--diversity"), dest="diversity", type="numeric", default=1,
					help="The diversity fraction. From 0 to 1. By default = 1.")

parser = add_option(parser, c("-r", "--replicates"), dest="replicates", type="numeric", default=10,
					help="The number of times to replicate the analysis. By default = 10")

parser = add_option(parser, c("-l", "--threshold"), dest="limit", type="numeric", default=5,
					help="The minimum number of times that a given node has a shift to be considered. By default = 5")

parser = add_option(parser, c("-p", "--prefix"), dest="prefix", type="character", default="NONE",
					help="The prefix for the output. By default, will add '_medusa_shifts.tsv' and '_medusa_summary.tsv' to the input tree, after deleting the extension.")

parser = add_option(parser, c("-c", "--cores"), dest="cores", type="numeric", default=1,
					help="The number of cores. By default = 1.")

parser = add_option(parser, c("-n", "--netDiversification"), dest="netDiversification", type="numeric", default=0.05,
					help="Initial conditions for net diversification. By default = 0.05.")

parser = add_option(parser, c("-e", "--relativeExtinction"), dest="relativeExtinction", type="numeric", default=0.5,
					help="Initial conditions for relative extinction By default = 0.5.")

parser = add_option(parser, c("-v", "--verbose"), dest="verbose", action="store_true", default=TRUE,
					help="If selected, will not print information to the console.")

args = parse_args(parser)

# Setting file names _______________________________________________________________________________
if(args$prefix == "NONE"){
	prefix = sub("\\.[^\\.]+$", "", args$tree)
}else{
	prefix = args$prefix
}

# Read files _______________________________________________________________________________________
if(args$verbose){cat("\n  Reading tree", args$tree, "\n")}
tree = read.tree(args$tree)
if(args$verbose){cat("  Diversity fraction", args$diversity, "\n")}
if(args$verbose){cat("  Replicates", args$replicates, "\n")}
if(args$verbose){cat("  Threshold", args$limit, "\n")}
if(args$verbose){cat("  Initial net diversification", args$netDiversification, "\n")}
if(args$verbose){cat("  Initial relative extinction", args$relativeExtinction, "\n")}
if(args$verbose){cat("  Results will be exported to", paste0(prefix, "_medusa_shifts.tsv"), "\n")}

# Starting analysis ________________________________________________________________________________
if(args$verbose){cat("  Starting analysis\n")}
medusa_shifts = data.frame()
for(rep in 1:args$replicates){
	if(args$verbose){cat("  Generating random tips")}
	richness = data.frame()
	for(i in 1:Ntip(tree)){
		cat("\r  Generating random tips ", round(i/Ntip(tree)*100), "%", sep="")
		r = 0
		while(r < 1){r = round(rnorm(1, 1/args$diversity, 0.5))}
		# r = round(rnorm(1, 1/args$diversity, 0.5))
		tmp = data.frame(tree$tip.label[i], ifelse(r<1, 1, r))
		richness = rbind(richness, tmp)
	}; rm(i, r, tmp)
	colnames(richness) = c("taxon", "n.taxa")
	if(args$verbose){cat("\n\n  Replicate ", rep, ": there are ", round(Ntip(tree)/args$diversity), " theoretical tips and ", sum(richness$n.taxa)," randomly generated tips\n\n", sep="")}
	res = medusa(tree, richness, model="bd", verbose=args$verbose,
				 init = c(r=args$netDiversification, epsilon=args$relativeExtinction),
				 ncores=args$cores)
	tmpr = res$summary
	tmpr$replicate = rep
	medusa_shifts = rbind(medusa_shifts, tmpr)
	if(rep == 1){
		write.table(tmpr, paste0(prefix, "_medusa_shifts.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	}else{
		write.table(tmpr, paste0(prefix, "_medusa_shifts.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
	}
}; rm(rep, tmpr)

# Summarising ______________________________________________________________________________________
if(args$verbose){cat("  Summarising analysis\n")}
nodes = as.data.frame(table(medusa_shifts$Shift.Node))
colnames(nodes) = c("nodes", "count")
nodes$pass = fifelse(nodes$count >= args$limit, TRUE, FALSE)

medusa_shifts$count = NA
medusa_shifts$pass = NA
medusa_shifts$age = NA
medusa_shifts$root = NA
root = max(nodeHeights(tree))
for(i in unique(medusa_shifts$Shift.Node)){
	medusa_shifts$count[which(medusa_shifts$Shift.Node==i)]=subset(nodes, nodes==i)$count
	medusa_shifts$pass[which(medusa_shifts$Shift.Node==i)]=subset(nodes, nodes==i)$pass
	medusa_shifts$age[which(medusa_shifts$Shift.Node==i)]=root-nodeheight(tree, i)
	medusa_shifts$root[which(medusa_shifts$Shift.Node==i)]=ifelse(nodeheight(tree, i) == 0, "root", "node")
}; rm(i)

heights = select(medusa_shifts, c("Shift.Node", "pass", "count", "age", "root"))
heights = heights[!duplicated(heights),]
heights = heights[order(-heights$count, -heights$age),]

# Exporting ________________________________________________________________________________________
if(args$verbose){cat("  Exporting final results file to", paste0(prefix, "_medusa_shifts.tsv"), "\n")}
write.table(medusa_shifts, paste0(prefix, "_medusa_shifts.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
if(args$verbose){cat("  Exporting summary to", paste0(prefix, "_medusa_summary.tsv"), "\n")}
write.table(heights, paste0(prefix, "_medusa_summary.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
if(args$verbose){cat("Done\n")}
