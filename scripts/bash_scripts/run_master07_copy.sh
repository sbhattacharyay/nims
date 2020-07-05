#!/bin/bash
#SBATCH --job-name=nims_shubs
#SBATCH --time=24:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=end
#SBATCH --mail-user=sbhatt15@jhu.edu
#
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------
 
ml matlab/R2019b-v2
ml # confirm modules used
matlab -nodisplay -nojvm -nosplash -nodesktop -r "cd('~/work/nims/scripts/'); master07_reverse;"
echo "matlab exit code: $?"