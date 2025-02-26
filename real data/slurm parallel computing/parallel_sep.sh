#!/bin/bash

#SBATCH -o Results_sep/Run_%a.Rout
#SBATCH --account=shenx
#SBATCH --array=1-50
#SBATCH --mail-user=pyrrhanikoss@outlook.com
#SBATCH --mail-type=END
#SBATCH --job-name=MSI_sep
#SBATCH --mem-per-cpu=1gb
#SBATCH --nodes=1
#SBATCH -t 10:00:00


module load R/4.4.0-openblas-rocky8

R CMD BATCH --vanilla main_sep.R  Results_sep/Run_${SLURM_ARRAY_TASK_ID}.Rout
