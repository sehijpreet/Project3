#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=E+L+G+P_orig
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skaur4@ufl.edu
#SBATCH --error=E+L+G+P_orig.err
#SBATCH --output=E+L+G+P_orig.out
#SBATCH --array=1-10  

module load R/4.2

Rscript code.R $SLURM_ARRAY_TASK_ID