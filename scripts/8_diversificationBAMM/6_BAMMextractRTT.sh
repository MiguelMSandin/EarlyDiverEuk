#!/bin/bash

# 'cd' to the root directory of a phylogenetic tree, where you have your dated tree, the subclades and the control files
# And now we will go to the 'bamm/diver' directory and extract independently every Rate Through Time plot and best shift configuration dates from the *event_data* files
cd "bamm/diver"
# The reason of running clade by clade is because some subclades are so large that it is required quite a lot of RAM memory, so we don't ask for 350GB of RAM all the time
# And the reason why these scripts are so convoluted (i.e.; this script calls '6.1_BAMMextractRTT_launcher.sh' and the latter calls '6.2_rttBAMMextract.R'),
# it is because by this way it is easier to deal with failures in the submitted jobs, launch some analyses for some clades and not others (some clades might finish within 2 hours, others might take up to 4 days), change headings regarding different cluster specifications, etc...

# In slurm clusters, you might want to call all the future jobs with a given prefix. For example:
# PREFIX="C11BE"

# This is the script we will call to run the bamm analyses
SCRIPT="PATH/TO/SCRIPT/6.1_BAMMextractRTT_launcher.sh"

cd Alveolata
mv ../*Alveolata* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Alveolata.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}alv$i --partition long --mem 240GB -t 4-0:00:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd ../Amoebozoa
mv ../*Amoebozoa* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Amoebozoa.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}amo$i --partition fast --mem 20GB -t 10:00:00 tmp.sh
done

cd ../Archaeplastida
mv ../*Archaeplastida* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Archaeplastida.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}arc$i --partition fast --mem 100GB -t 1-0:00:00 tmp.sh
done

cd ../Cryptista
mv ../*Cryptista* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Cryptista.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}cry$i --partition fast --mem 10GB -t 10:00:00 tmp.sh
done

cd ../Discoba
mv ../*Discoba* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Discoba.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}dis$i --partition fast --mem 80GB -t 1-0:00:00 tmp.sh
done

cd ../Haptista
mv ../*Haptista* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Haptista.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}hap$i --partition fast --mem 10GB -t 10:00:00 tmp.sh
done

cd ../Holozoa
mv ../*Holozoa* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Holozoa.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}hol$i --partition bigmem --mem 350GB -t 4-0:00:00 tmp.sh
done

cd ../Metamonada
mv ../*Metamonada* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Metamonada.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}met$i --partition fast --mem 10GB -t 10:00:00 tmp.sh
done

cd ../Nucletmycea
mv ../*Nucletmycea* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Nucletmycea.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}nuc$i --partition bigmem --mem 300GB -t 4-0:00:00 tmp.sh
done

cd ../Rhizaria
mv ../*Rhizaria* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Rhizaria.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}rhi$i --partition fast --mem 120GB -t 1-0:00:00 tmp.sh
done

cd ../Stramenopila
mv ../*Stramenopila* .
i=0
for FILE in $(ls *event_data*); do
	((i=i+1))
	TREE="../../../clades/clade_Stramenopila.tre"
	cp $SCRIPT tmp.sh
	sed -i "s|INSERT_TREE_HERE|$TREE|g" tmp.sh
	sed -i "s|INSERT_EDATA_HERE|$FILE|g" tmp.sh
	sbatch --job-name ${PREFIX}str$i --partition fast --mem 120GB -t 1-0:00:00 tmp.sh
done

cd ../
tree -h
#! bin/bash
