#!/bin/bash

#SBATCH -o Results_CE/Run_%a.Rout
#SBATCH --array=1-25
#SBATCH --mail-user=pyrrhanikoss@outlook.com
#SBATCH --mail-type=END
#SBATCH --job-name=MSI_cross_entropy
#SBATCH --mem-per-cpu=1gb
#SBATCH --nodes=1
#SBATCH -t 24:00:00

module load R

R CMD BATCH --vanilla main_CE.R  Results_CE/Run_${SLURM_ARRAY_TASK_ID}.Rout
