#!/bin/bash

# 'cd' to the 'clads' directory created with the previous script
# The table with the minimum or maximum diversity fractions for every clade
FRACTIONS="fractions_max.tsv"

# In slurm clusters, you might want to call all the future jobs starting with a given prefix. For example:
# PREFIX="C11C"

# This is the script we will call to run the CLaDS analyses
SCRIPT="PATH/TO/SCRIPT/1.1_CLaDS_launcher.sh"

# Alveolata
FRACTION=$(grep Alveolata $FRACTIONS | cut -f 3)
TREE=$(grep Alveolata $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Amoebozoa
FRACTION=$(grep Amoebozoa $FRACTIONS | cut -f 3)
TREE=$(grep Amoebozoa $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Archaeplastida
FRACTION=$(grep Archaeplastida $FRACTIONS | cut -f 3)
TREE=$(grep Archaeplastida $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Cryptista
FRACTION=$(grep Cryptista $FRACTIONS | cut -f 3)
TREE=$(grep Cryptista $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Discoba
FRACTION=$(grep Discoba $FRACTIONS | cut -f 3)
TREE=$(grep Discoba $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Haptista
FRACTION=$(grep Haptista $FRACTIONS | cut -f 3)
TREE=$(grep Haptista $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Holozoa
FRACTION=$(grep Holozoa $FRACTIONS | cut -f 3)
TREE=$(grep Holozoa $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Metamonada
FRACTION=$(grep Metamonada $FRACTIONS | cut -f 3)
TREE=$(grep Metamonada $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Nucletmycea
FRACTION=$(grep Nucletmycea $FRACTIONS | cut -f 3)
TREE=$(grep Nucletmycea $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Rhizaria
FRACTION=$(grep Rhizaria $FRACTIONS | cut -f 3)
TREE=$(grep Rhizaria $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh

# Stramenopila
FRACTION=$(grep Stramenopila $FRACTIONS | cut -f 3)
TREE=$(grep Stramenopila $FRACTIONS | cut -f 1)
FRAC=$(echo $FRACTION | awk '{printf("%.2f",$1)}' | sed 's/0\.//g')
cp $SCRIPT tmp.sh
sed -i "s/INSERT_TREE_HERE/$TREE/g" tmp.sh
sed -i "s/INSERT_FRAC_HERE/$FRACTION/g" tmp.sh
# In slurm cluster something like this will do:
# sbatch --mem 80GB -t 20-0:00 --job-name=$PREFIX"alv"$FRAC tmp.sh
# Or instead run it locally as:
bash tmp.sh
