#!/bin/bash

#SBATCH --job-name=func1
#SBATCH --time=0-00:10:00
#SBATCH --output=func_mqtls_birth_450k_cent5.out
# #SBATCH --mem=20GB
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl


Rscript 3.functional_analysis1_mqtls.R 'within_centile' &
Rscript 3.functional_analysis1_mqtls.R 'across_centiles' &
wait

exit 0
