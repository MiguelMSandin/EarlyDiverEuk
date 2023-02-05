#! bin/bash

# 'cd' to the root directory of a phylogenetic tree, where you have your dated tree, the subclades and the control files

# We create a subdirectory in the 'bamm' directory to store the output analyses from BAMM
[ ! -d "bamm/diver" ] && mkdir -p "bamm/diver"

# In slurm clusters, you might want to call all the future jobs with a given prefix. For example:
# PREFIX="C11Bt"

# This is the script we will call to run the bamm analyses
SCRIPT="PATH/TO/SCRIPT/2.1_BAMMestimateBestShifts_launcher.sh"

# And now we will go to the 'bamm/control_test' directory and run independently every control file
cd "bamm/control_test"

cd "Alveolata"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}alv$i --mem 20GB -t 5-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Amoebozoa"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}amo$i --mem 4GB -t 2-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Archaeplastida"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}arc$i --mem 6GB -t 4-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Cryptista"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}cry$i --mem 2GB -t 2-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Discoba"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}dis$i --mem 6GB -t 4-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Haptista"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}hap$i --mem 2GB -t 2-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Holozoa"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}hol$i --mem 20GB -t 10-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Metamonada"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}met$i --mem 4GB -t 2-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Nucletmycea"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}nuc$i --mem 20GB -t 10-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Rhizaria"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}rhi$i --mem 20GB -t 10-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../Stramenopila"
i=0
for FILE in $(ls *ctl); do
	((i=i+1))
	cp "$SCRIPT" "tmp.sh"
	sed -i "s/INSERT_FILE_NAME_HERE/$FILE/g" tmp.sh
	# In slurm cluster something like this will do:
	# sbatch --job-name ${PREFIX}str$i --mem 20GB -t 10-0:00 tmp.sh
	# Or instead run it locally as:
	bash tmp.sh
done

cd "../"
