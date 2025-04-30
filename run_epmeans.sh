#!/bin/bash

#SBATCH --job-name=epmeans2
#SBATCH --time=3-00:00:00
#SBATCH --output=epmeans_birth_450k_cent5.out
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=END,FAIL


Rscript 1.fit_epmeans.R

exit 0
