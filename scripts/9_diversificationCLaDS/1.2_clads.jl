
using ArgParse;
using PANDA;
using JLD2;

function parse_commandline()
	s = ArgParseSettings()
	@add_arg_table s begin
		"--tree", "-t"
			help = "A time calibrated phylogenetic tree, in newick format, to be computed the ClaDS function."
			required = true
		"--fraction", "-f"
			help = "The fraction of the diversity sampled. Deafult = 1"
			#arg_type = Float64
			default = 1
		"--outputRdata", "-r"
			help = "The output file saved to Rdata."
			required = true
# 		"--outputJulia", "-j"
# 			help = "The output file saved to Julia."
# 			required = true
#		"--verbose", "-v""
#			help = "an option without argument"
#			action = :store_true
	end
	return parse_args(s)
end

function run_clads()
	parsed_args = parse_commandline();
	tree = load_tree(parsed_args["tree"]);
	fraction = parse(Float64, parsed_args["fraction"]);
	clads = infer_ClaDS(tree, f = fraction, print_state = 100);
# 	@save parsed_args["outputJulia"] clads;
# 	plot_CladsOutput(clads, method = "tree");
# 	plot_CladsOutput(clads, method = "DTT");
# 	plot_CladsOutput(clads, method = "RTT");
	save_ClaDS_in_R(clads, parsed_args["outputRdata"]);
end

run_clads()
