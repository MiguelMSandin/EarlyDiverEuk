#!/bin/bash

# mkdir medusa
# cd medusa
# cp ../clads/fractions.tsv .
# cp ../clades/*tre .

# bash /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh

# The fractions file must have 2 columns (tree\tfraction) and in 0 to 1 format!
FRACTIONS="fractions.tsv"
PREFIX="Dfc22M"

# Alveolata
i=0
for FRACTION in $(grep Alveolata $FRACTIONS | cut -f 2); do
	TREE="clade_Alveolata.tre"
	((i=i+1))
	JOBNAME=$PREFIX"alv"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Alveolata/g" tmp.sh

	sbatch --partition long --mem 6GB -t 6-0:00 --job-name=$JOBNAME tmp.sh
done

# Amoebozoa
i=0
for FRACTION in $(grep Amoebozoa $FRACTIONS | cut -f 2); do
	TREE="clade_Amoebozoa.tre"
	((i=i+1))
	JOBNAME=$PREFIX"amo"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Amoebozoa/g" tmp.sh

	sbatch --partition fast --mem 2GB --partition fast -t 0-4:00:00 --job-name=$JOBNAME tmp.sh
done

# Archaeplastida
i=0
for FRACTION in $(grep Archaeplastida $FRACTIONS | cut -f 2); do
	TREE="clade_Archaeplastida.tre"
	((i=i+1))
	JOBNAME=$PREFIX"arc"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Archaeplastida/g" tmp.sh

	sbatch --partition fast --mem 3GB -t 1-0:00 --job-name=$JOBNAME tmp.sh
done

# Cryptista
i=0
for FRACTION in $(grep Cryptista $FRACTIONS | cut -f 2); do
	TREE="clade_Cryptista.tre"
	((i=i+1))
	JOBNAME=$PREFIX"cry"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Cryptista/g" tmp.sh

	sbatch --partition fast --mem 2GB --partition fast -t 1:00:00 --job-name=$JOBNAME tmp.sh
done

# Discoba
i=0
for FRACTION in $(grep Discoba $FRACTIONS | cut -f 2); do
	TREE="clade_Discoba.tre"
	((i=i+1))
	JOBNAME=$PREFIX"dis"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Discoba/g" tmp.sh

	sbatch --partition fast --mem 2GB --partition fast -t 1-0:00 --job-name=$JOBNAME tmp.sh
done

# Haptista
i=0
for FRACTION in $(grep Haptista $FRACTIONS | cut -f 2); do
	TREE="clade_Haptista.tre"
	((i=i+1))
	JOBNAME=$PREFIX"hap"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Haptista/g" tmp.sh

	sbatch --partition fast --mem 2GB --partition fast -t 0-2:00:00 --job-name=$JOBNAME tmp.sh
done

# Holozoa
i=0
for FRACTION in $(grep Holozoa $FRACTIONS | cut -f 2); do
	TREE="clade_Holozoa.tre"
	((i=i+1))
	JOBNAME=$PREFIX"hol"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Holozoa/g" tmp.sh

	sbatch --partition long --mem 6GB -t 10-0:00 --job-name=$JOBNAME tmp.sh
done

# Metamonada
i=0
for FRACTION in $(grep Metamonada $FRACTIONS | cut -f 2); do
	TREE="clade_Metamonada.tre"
	((i=i+1))
	JOBNAME=$PREFIX"met"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Metamonada/g" tmp.sh

	sbatch --partition fast --mem 2GB --partition fast -t 01:00:00 --job-name=$JOBNAME tmp.sh
done

# Nucletmycea
i=0
for FRACTION in $(grep Nucletmycea $FRACTIONS | cut -f 2); do
	TREE="clade_Nucletmycea.tre"
	((i=i+1))
	JOBNAME=$PREFIX"nuc"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Nucletmycea/g" tmp.sh

	sbatch --partition long --mem 4GB -t 6-0:00 --job-name=$JOBNAME tmp.sh
done

# Rhizaria
i=0
for FRACTION in $(grep Rhizaria $FRACTIONS | cut -f 2); do
	TREE="clade_Rhizaria.tre"
	((i=i+1))
	JOBNAME=$PREFIX"rhi"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Rhizaria/g" tmp.sh

	sbatch --partition fast --mem 3GB -t 1-0:00 --job-name=$JOBNAME tmp.sh
done

# Stramenopila
i=0
for FRACTION in $(grep Stramenopila $FRACTIONS | cut -f 2); do
	TREE="clade_Stramenopila.tre"
	((i=i+1))
	JOBNAME=$PREFIX"str"$i

	cp /shared/projects/radecoevo/pacbio/scripts/job_geiger_launcher.sh tmp.sh
	sed -i "s/INSERT_TREE/$TREE/g" tmp.sh
	sed -i "s/INSERT_FRACTION/$FRACTION/g" tmp.sh
	sed -i "s/INSERT_CLADE/Stramenopila/g" tmp.sh

	sbatch --partition fast --mem 3GB -t 1-0:00 --job-name=$JOBNAME tmp.sh
done

rm -f tmp.sh
