#!/usr/bin/env python3

import argparse
from Bio import Phylo
import re

parser = argparse.ArgumentParser(description="Given a control file and a directory of trees, export how many tips contain each calibration per tree.")

# Add the arguments to the parser
requiredArgs = parser.add_argument_group('required arguments')

requiredArgs.add_argument("-c", "--controlFile", dest="controlFile", required=True,
					help="A control file for treePL.")

requiredArgs.add_argument("-t", "--trees", dest="trees", required=True, nargs="+",
						  help="A list of tree files in newick format.")

parser.add_argument("-o", "--out", dest="out", required=False, default=None,
					help="The name of the output tsv file. By default, will add '_summary.tsv' to the control file.")

parser.add_argument("-v", "--verbose", dest="verbose", required=False, action="store_false",
					help="If selected, will not print information to the console.")

args = parser.parse_args()

# Reading files ------------------------------------------------------------------------------------
if args.verbose:
	print("  Reading control file:", args.controlFile)
taxa = {}
taxaList = list()
agesMin = {}
agesMax = {}
for line in open(args.controlFile):
	if line.startswith("mrca"):
		line = line.strip(" ").split()
		name = str(line[2])
		taxa[name] = list()
		for taxon in range(3, len(line)):
			taxa[name].append(line[taxon])
			taxaList.append(line[taxon])
	elif line.startswith("min"):
		line = line.strip(" ").split()
		agesMin[str(line[2])] = str(line[3])
	elif line.startswith("max"):
		line = line.strip(" ").split()
		agesMax[str(line[2])] = str(line[3])

if args.out is None:
	out = re.sub("\\.[^\\.]+$", "_summary.tsv", args.controlFile)
else:
	out = args.out

# Print a quick information on the stuff to be analyzed --------------------------------------------
if args.verbose:
	print("    There are", len(taxa), "calibration nodes")

# Check number of tips for every calibratated node -------------------------------------------------
if args.verbose:
	print("  Analyzing tree", end="")
	length = len(args.trees)
	i = 0
with open(out, "w") as outfile:
	tmp = "\t".join(taxa.keys())
	print("tree\t" + tmp, file=outfile)
	for tree in args.trees:
		tmp = list()
		if args.verbose:
			i += 1
			print("\r  Analyzing tree ", i, "/", length, sep="", end="")
		T = Phylo.read(tree, "newick")
		for tips in taxa.values():
			clade = T.common_ancestor(tips)
			count = clade.count_terminals()
			tmp.append(str(count))
		tmp = "\t".join(tmp)
		print(tree + "\t" + tmp, file=outfile)

if args.verbose:
	print("\nDone")
