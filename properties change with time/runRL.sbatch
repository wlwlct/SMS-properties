#!/bin/bash

#SBATCH -J $jobname
#SBATCH -p interactive
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu=32gb
#SBATCH --cpus-per-task=1
#SBATCH -e error.err
#SBATCH -o output.out
module load matlab
matlab -singleCompThread -nodisplay -nosplash -r propertieschangewithtimeRmLI
##-nojvm
