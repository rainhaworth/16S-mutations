#!/bin/bash
#SBATCH --output=~/16S-mutations/results/run.o
#SBATCH --time=13:00:00
#SBATCH --mem=32gb

source ~/.bashrc; conda activate tfpy39

srun python gapped_kmers.py -k 8 -n 600000 -t 0.95 --save
