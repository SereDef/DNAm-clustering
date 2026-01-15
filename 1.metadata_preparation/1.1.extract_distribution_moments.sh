#!/bin/bash

#SBATCH --job-name=1.1.moments
#SBATCH --time=1-00:00:00
#SBATCH --output=1.1.log
# #SBATCH --mem=50GB
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl

# Run on the mega matrix (all cohorts all arrays) 
# before and after batch correction (sequentially)
Rscript 1.1.extract_distribution_moments.R "clean"
Rscript 1.1.extract_distribution_moments.R "batch_corrected"

# Run the two arrays (sequentially)
# Rscript 1.1.extract_distribution_moments.R "450K"
# Rscript 1.1.extract_distribution_moments.R "EPIC"


