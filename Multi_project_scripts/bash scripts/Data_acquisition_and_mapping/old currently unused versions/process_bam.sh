#!/bin/bash

#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="_%x_%j.out"
#SBATCH --mail-type=all

#For SBATCH commands above, the x gives the job name and the j gives the job ID.
#module load bamtools
#module load deeptools

#for f in *.bam;
#       do
#       samtools index $f #generates bam.bai for every bam file in the folder
#       outFile=$(basename "$f" .bam | awk '{print ""$1".bw"}') #makes a file name for output. Takes filename - .bam, then adds .bw to that string.
#        bamCoverage --bam $f \
#               --outFileName $outFile \
#               --binSize 5 \
#               --normalizeUsingRPKM ##makes bigwig
#done


