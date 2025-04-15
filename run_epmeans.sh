#!/bin/bash

#SBATCH --job-name=epmeans
#SBATCH --time=4-12:00:00
#SBATCH --output=birth_epmeans_cent5
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL


Rscript 1.fit_epmeans.R

exit 0
