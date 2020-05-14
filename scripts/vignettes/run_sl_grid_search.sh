#!/bin/bash

#!/bin/bash
#SBATCH --job-name=rscript
#SBATCH --time=0:0:10
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
#---------------------------------------------------------------------
# SLURM job script to run serial R
#---------------------------------------------------------------------
 
ml R
ml # confirm modules used
Rscript sl_grid_search.R
