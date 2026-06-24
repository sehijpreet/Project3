#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=ssBLUP
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skaur4@ufl.edu
#SBATCH --error=ssBLUP.err
#SBATCH --output=ssBLUP.out
#SBATCH --array=1-10  

module load R/4.2

Rscript code1.R $SLURM_ARRAY_TASK_ID