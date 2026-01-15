#!/bin/bash

#SBATCH --job-name=2.3.epimeans
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --error=2.3.log
#SBATCH --output=2.3.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl

Rscript 2.3.fit_epmeans_stage2.R
