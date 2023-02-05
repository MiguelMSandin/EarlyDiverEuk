#! bin/bash

# Run 10 independent phylogenetic analyses with 10 bootstrap each --------------------------------
FASTA="pb_pr2_reads10c_align_trim05.fasta"
CONSTRAINT="eukProt_v9_forConstraint.tre"
BS="10"
THREADS="16"

# In slurm clusters you might want to change this loop for '#SBATCH --array=1-10'
# About 30GB of RAM memory and 6 days will be sufficient
for SLURM_ARRAY_TASK_ID in $(seq 10); do
	OUTPUT=${FASTA/.fasta/_raxml-GTRgamma-D$SLURM_ARRAY_TASK_ID}
	raxmlHPC-PTHREADS-AVX2 -g $CONSTRAINT -T $THREADS -m GTRGAMMA -p $RANDOM -D -N $BS -n $OUTPUT -s $FASTA --no-seq-check
done

# Now search for the best scoring tree and annotate the 100 total bootstrap in the tree ------------
BS_ALL="RAxML_bootstrap.${FASTA/.fasta/_raxml-GTRgamma.tre}"
OUT=${FASTA/.fasta/_raxml-GTRgamma}

# Concatenating all bootstraps into one single file
for BS in $(ls | grep "RAxML_bootstrap\.${FASTA/.fasta/}"); do
    cat $BS >> $BS_ALL
done

# Searching best tree
for INFO in $(ls | grep "RAxML_info\.${FASTA/.fasta/}"); do
    like=$(grep "of best tree " $INFO | sed -e 's/.* //g')
    echo "$INFO: $like" >> "$OUT.likelihoods"
    echo "$INFO: $like"
done

MIN=$(cut -f2 -d ":" "${OUT/.tre/.likelihoods}" | sort -n | tail -1)
BEST_TREE=$(grep "$MIN" "${OUT/.tre/.likelihoods}" | sed -e 's/: .*//g')
BEST_TREE=${BEST_TREE/RAxML_info/RAxML_bestTree}

# Drawing bipartitions
module load raxml/8.2.12
raxmlHPC-PTHREADS-SSE3 -m GTRCAT -p $(date +%s) -f b -t $BEST_TREE -z $BS_ALL -n $OUT

echo ""
echo "Best tree found in: $BEST_TREE"
