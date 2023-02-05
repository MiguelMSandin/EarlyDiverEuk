#! bin/bash

# Download in-house scripts ------------------------------------------------------------------------
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/fasta/fastaRevCom.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/treePruneOutliers.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/treeRootOutgroup.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/treeColourBranches.py
wget https://raw.githubusercontent.com/MiguelMSandin/random/main/phylogenetics/treeRemoveBranchLengths.py

chmod +x *py
mv *py ../

# Align forward dataset ----------------------------------------------------------------------------
# In slurm cluster you might ask for 100GB of RAM memory and about 5 days to be in the safe side
THREADS="16"
mafft --quiet --thread $THREADS "pb_pr2_reads10c.fasta" > "pb_pr2_reads10c_align.fasta"

# Trim
# In slurm cluster you might ask for 100GB of RAM memory and 1 day might be enough
trimal -in "pb_pr2_reads10c_align.fasta" -out "pb_pr2_reads10c_align_trim05.fasta" -gt "0.05"
gzip "pb_pr2_reads10c_align.fasta"

# Align reverse dataset ----------------------------------------------------------------------------
# Reverse
fastaRevCom.py -f "pb_pr2_reads10c.fasta" -o "pb_pr2_reads10c_rev.fasta"

# Align
# In slurm cluster you might ask for 100GB of RAM memory and about 5 days to be in the safe side
mafft --quiet --thread $THREADS "pb_pr2_reads10c_rev.fasta" > "pb_pr2_reads10c_rev_align.fasta"

# Trim
# In slurm cluster you might ask for 100GB of RAM memory and 1 day might be enough
trimal -in "pb_pr2_reads10c_rev_align.fasta" -out "pb_pr2_reads10c_rev_align_trim05.fasta" -gt "0.05"
gzip "pb_pr2_reads10c_rev_align.fasta"

# Reverse again
fastaRevCom.py -f "pb_pr2_reads10c_rev_align_trim05.fasta" -o "pb_pr2_reads10c_r2_align_trim05.fasta"
gzip "pb_pr2_reads10c_rev_align_trim05.fasta"

# And finally organise the datasets into different folders -----------------------------------------
mkdir step1f step1r
mv *rev* step1r
mv pb_pr2* step1f
