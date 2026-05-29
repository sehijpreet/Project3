#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=blues_pheno
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skaur4@ufl.edu
#SBATCH --error=blues_Pheno.err
#SBATCH --output=blues_Pheno.out
#SBATCH --array=1-5  

module load R/4.2

Rscript code_p.R $SLURM_ARRAY_TASK_ID