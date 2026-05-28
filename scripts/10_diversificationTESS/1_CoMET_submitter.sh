#!/bin/bash

# mkdir comet
# cd comet
# cp ../clades/clade_* .

# bash /shared/projects/radecoevo/pacbio/scripts/comet_submitter.sh

# The fractions file must have 2 columns (tree\tfraction) and in 0 to 1 format!
FRACTIONS="fractions.tsv"
PREFIX="2Ar621C"
LAUNCHER="/shared/projects/radecoevo/pacbio/scripts/comet_launcher.sh"

# Alveolata
i=0
for FRACTION in $(grep Alveolata $FRACTIONS | cut -f 2); do
	TREE="clade_Alveolata.tre"
	((i=i+1))
	JOBNAME=$PREFIX"alv"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 20GB -t 0-6:00:00 --job-name=$JOBNAME tmp.sh
done

# Amoebozoa
i=0
for FRACTION in $(grep Amoebozoa $FRACTIONS | cut -f 2); do
	TREE="clade_Amoebozoa.tre"
	((i=i+1))
	JOBNAME=$PREFIX"amo"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 1GB -t 0-1:00:00 --job-name=$JOBNAME tmp.sh
done

# Archaeplastida
i=0
for FRACTION in $(grep Archaeplastida $FRACTIONS | cut -f 2); do
	TREE="clade_Archaeplastida.tre"
	((i=i+1))
	JOBNAME=$PREFIX"arc"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 10GB -t 0-2:00:00 --job-name=$JOBNAME tmp.sh
done

# Cryptista
i=0
for FRACTION in $(grep Cryptista $FRACTIONS | cut -f 2); do
	TREE="clade_Cryptista.tre"
	((i=i+1))
	JOBNAME=$PREFIX"cry"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 1GB -t 00:30:00 --job-name=$JOBNAME tmp.sh
done

# Discoba
i=0
for FRACTION in $(grep Discoba $FRACTIONS | cut -f 2); do
	TREE="clade_Discoba.tre"
	((i=i+1))
	JOBNAME=$PREFIX"dis"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 3GB -t 2:00:00 --job-name=$JOBNAME tmp.sh
done

# Haptista
i=0
for FRACTION in $(grep Haptista $FRACTIONS | cut -f 2); do
	TREE="clade_Haptista.tre"
	((i=i+1))
	JOBNAME=$PREFIX"hap"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 1GB -t 0-1:00:00 --job-name=$JOBNAME tmp.sh
done

# Holozoa
i=0
for FRACTION in $(grep Holozoa $FRACTIONS | cut -f 2); do
	TREE="clade_Holozoa.tre"
	((i=i+1))
	JOBNAME=$PREFIX"hol"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 30GB -t 0-10:00:00 --job-name=$JOBNAME tmp.sh
done

# Metamonada
i=0
for FRACTION in $(grep Metamonada $FRACTIONS | cut -f 2); do
	TREE="clade_Metamonada.tre"
	((i=i+1))
	JOBNAME=$PREFIX"met"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 1GB -t 00:30:00 --job-name=$JOBNAME tmp.sh
done

# Nucletmycea
i=0
for FRACTION in $(grep Nucletmycea $FRACTIONS | cut -f 2); do
	TREE="clade_Nucletmycea.tre"
	((i=i+1))
	JOBNAME=$PREFIX"nuc"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 20GB -t 0-6:00:00 --job-name=$JOBNAME tmp.sh
done

# Rhizaria
i=0
for FRACTION in $(grep Rhizaria $FRACTIONS | cut -f 2); do
	TREE="clade_Rhizaria.tre"
	((i=i+1))
	JOBNAME=$PREFIX"rhi"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 10GB -t 0-2:00:00 --job-name=$JOBNAME tmp.sh
done

# Stramenopila
i=0
for FRACTION in $(grep Stramenopila $FRACTIONS | cut -f 2); do
	TREE="clade_Stramenopila.tre"
	((i=i+1))
	JOBNAME=$PREFIX"str"$i
	cp ${LAUNCHER} tmp.sh
	sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
	sbatch --partition fast --mem 10GB -t 0-2:00:00 --job-name=$JOBNAME tmp.sh
done

rm -f tmp.sh

