#!/bin/bash

#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="import_%u_%j.out"
#SBATCH --mail-type=all

cat skew.tr3.bed | awk 'BEGIN {srand()} !/^$/ {if (rand() <= .01) print $0}' > subset.bed
./bedGraphToBigWig subset.bed chromSizes.txt subset.bw
