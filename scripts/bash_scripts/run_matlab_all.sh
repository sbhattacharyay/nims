#!/bin/bash
#SBATCH --job-name=all_matlab_nims
#SBATCH --time=72:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
### Using more tasks because default memory is ~5GB per core
### 'shared' will share the node with other users
### 'parallel' use entire node (24,28,48, depends on node type)
### Try the script with --ntasks-per-node=1 and see what happens
#
#---------------------------------------------------------------------
# SLURM job script to run serial MATLAB
#---------------------------------------------------------------------
 
ml matlab/R2019b-v2
ml # confirm modules used
matlab -nodisplay -nojvm -nosplash -nodesktop -r "cd('~/work/SB_bims/scripts/'); master02_all_motion_feature_extraction;"
echo "matlab exit code: $?"