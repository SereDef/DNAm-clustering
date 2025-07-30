#!/bin/bash

#SBATCH --job-name=0.3.combat
#SBATCH --time=1-00:00:00
#SBATCH --output=0.3.log
# #SBATCH --mem=50GB
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl


Rscript 0.3.filter_and_batch_correction.R

exit 0
