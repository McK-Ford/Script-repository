#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --mail-type=all

###Above tells SLURM this is a slurm script, that I want it to allocate 12G memory for, limit to running 12 hours, and email me when it starts/finishes/fails.

module load deeptools
for f in *.bam;
	do
	outFile=$(basename "$f" .bam | awk '{print $1".bw"}') #makes filename by taking bam off filename then adding bw to it.
	bamCoverage --bam $f \
		--outFileName $outFile \
		--binSize 10 #makes a bw file with 10bp bins
done
