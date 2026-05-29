#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=blues_res
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skaur4@ufl.edu
#SBATCH --error=blues_res.err
#SBATCH --output=blues_res.out
#SBATCH --array=1-5  

module load R/4.2

Rscript code_r.R $SLURM_ARRAY_TASK_ID