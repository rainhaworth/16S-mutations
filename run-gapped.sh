#!/bin/bash
#SBATCH --output=/nfshomes/rhaworth/16S-mutations/results/run-14-nopart.o
#SBATCH --time=12:00:00
#SBATCH --mem=32gb

source ~/.bashrc; conda activate /fs/nexus-scratch/rhaworth/env/tfpt

srun python gapped_kmers.py -k 14 -n 110000 -t 0.9 --save
