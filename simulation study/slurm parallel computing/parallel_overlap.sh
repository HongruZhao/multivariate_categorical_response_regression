#!/bin/bash

#SBATCH -o Results_overlap/Run_%a.Rout
#SBATCH --array=1-50
#SBATCH --mail-user=pyrrhanikoss@outlook.com
#SBATCH --mail-type=END
#SBATCH --job-name=MSI_overlap
#SBATCH --mem-per-cpu=1gb
#SBATCH --nodes=1
#SBATCH -t 24:00:00

module load R

R CMD BATCH --vanilla main_overlap.R  Results_overlap/Run_${SLURM_ARRAY_TASK_ID}.Rout
