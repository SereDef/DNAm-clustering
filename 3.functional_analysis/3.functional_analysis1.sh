#!/bin/bash

#SBATCH --job-name=3.func1
#SBATCH --array=1-2 # within and across range centiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64   # NOTE R limit parallel processes = 124
#SBATCH --time=3-00:00:00
#SBATCH --error=3.logs/3_log%a
#SBATCH --output=3.logs/3_log%a
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl

Rscript 3.functional_analysis1_mqtls.R