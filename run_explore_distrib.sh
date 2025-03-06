#!/bin/bash

#SBATCH --job-name=expdist
#SBATCH --output=bycentiles5
#SBATCH --mem=5GB
#SBATCH --mail-type=END,FAIL


Rscript MPSR_explore_distrib.R

exit 0
