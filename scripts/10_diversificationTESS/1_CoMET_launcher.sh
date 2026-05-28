#!/bin/bash
#
#SBATCH --output=slurm.%x.%N.%j.out
#SBATCH --error=slurm.%x.%N.%j.err
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user EMAILADDRESS
#SBATCH --cpus-per-task 1
#SBATCH --exclude=cpu-node-057

module load r/4.5.1

TREE="INSERT_TREE_HERE"
FRACTION="INSERT_FRAC_HERE"
FRAC=$(echo ${FRACTION} | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
CLADE=$(echo $TREE | sed "s/.*clade_//g" | sed "s/\.tre//g")
DIR="${CLADE}/d${FRAC}"
PREFIX="${CLADE}/${CLADE}_d${FRAC}"

echo "CoMET.R -t ${TREE} -f ${FRACTION} -d ${DIR} -p ${PREFIX}"

Rscript /shared/projects/radecoevo/pacbio/scripts/CoMET.R -t ${TREE} -f ${FRACTION} -d ${DIR} -p ${PREFIX}
