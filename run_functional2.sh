#!/bin/bash

#SBATCH --job-name=func2
#SBATCH --time=0-10:00:00
#SBATCH --output=func_ewascat_birth_450k_cent5.out
#SBATCH --mem=50GB
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl


Rscript 4.functional_analysis2_EWAScatalog.R 'within_centile' &
Rscript 4.functional_analysis2_EWAScatalog.R 'across_centiles' &
wait

exit 0
