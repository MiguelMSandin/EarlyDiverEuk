#! bin/bash

FASTA="pb_pr2_rev_align_trim05.fasta"
# Now this is starting to be confussing, so bear in mind what tree you are using for the constraint, and keep clear file names
CONSTRAINT=""
PREFIX=""

MEM="40G"
THREADS="16"

# In slurm clusters you might want to change this loop for '#SBATCH --array=1-2'
# About 40GB of RAM memory and 4 days will be sufficient to ensure completion
for SLURM_ARRAY_TASK_ID in $(seq 2); do
	OUTPUT=$PREFIX"_"$SLURM_ARRAY_TASK_ID
	iqtree -s $FASTA --seqtype DNA --prefix $OUTPUT -g $CONSTRAINT -m GTR+I+G --seed $RANDOM --mem $MEM -T $THREADS -fast
done

