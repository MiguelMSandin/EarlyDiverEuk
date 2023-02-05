#! bin/bash

# Run 10 independent phylogenetic analyses -------------------------------------------------------
FILE="pb_pr2_reads10c_r2_align_trim05.fasta"
CONSTRAINT="eukProt_v9_forConstraint.tre"
THREADS="16"

# In slurm clusters you might want to change this loop for '#SBATCH --array=1-10'
# About 30GB of RAM memory and 6 days will be sufficient
for SLURM_ARRAY_TASK_ID in $(seq 10); do
	PREFIX="step1r_RAng$SLURM_ARRAY_TASK_ID"
	raxml-ng --parse --msa $FILE --model GTR+G --prefix $PREFIX
	raxml-ng --msa $PREFIX.raxml.rba --tree-constraint $CONSTRAINT --tree pars{1} --prefix $PREFIX --seed $RANDOM --threads $THREADS
done

# Find the 4 best trees ----------------------------------------------------------------------------
grep "Final LogLikelihood" *log > "likelihoods_reverse.tsv"
sed -i "s/:Final LogLikelihood: /\t/g" "likelihoods_reverse.tsv"

for LL in $(cut -f 2 "likelihoods_reverse.tsv" | sort | head -n 4); do
	LLL=${LL/-/}
	BEST_TREE=$(grep "$LLL" "likelihoods_reverse.tsv" | cut -f 1)
	echo "$BEST_TREE"
	echo "$BEST_TREE" >> "best_trees_reverse.list"
done
