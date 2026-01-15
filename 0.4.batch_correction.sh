#!/bin/bash

#SBATCH --job-name=0.4.combat
#SBATCH --time=1-00:00:00
#SBATCH --output=0.4.log
# #SBATCH --mem=50GB
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl

# Run the on the mega matrix (all cohorts all arrays)
Rscript 0.4.batch_correction.R

# Run the two arrays (sequentially)
# Rscript 0.3.filter_and_batch_correction.R "EPIC" "GENR"
# Rscript 0.3.filter_and_batch_correction.R "450K" "GENR"
# Rscript 0.3.filter_and_batch_correction.R "450K" "ALSP"