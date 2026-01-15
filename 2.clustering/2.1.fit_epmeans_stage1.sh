#!/bin/bash

#SBATCH --job-name=2.1.epimeans
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:00
#SBATCH --error=2.1.logs/2.1_log%a
#SBATCH --output=2.1.logs/2.1_log%a
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=s.defina@erasmusmc.nl

Rscript 2.1.fit_epmeans_stage1.R
