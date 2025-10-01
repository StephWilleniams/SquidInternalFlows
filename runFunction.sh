#!/bin/bash

#SBATCH --job-name=squiderific
#SBATCH --time=0-06:00:00 
#SBATCH --array=1-3
#SBATCH --ntasks=40
#SBATCH --partition=short
#SBATCH --output=logs/matlab_job_%a.out
#SBATCH --error=logs/matlab_job_%a.err

# Set up MATLAB environment
module load matlab/r2024a

# Submit each array job
matlab -nodisplay -nodesktop -nosplash -r "CalculateTrajectory($SLURM_ARRAY_TASK_ID)"
