#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="Acquire_Fastqs_%j.out"
#SBATCH --mail-type=all

module load sratoolkit
awk '{print $0=$1}' Refsheet.txt > sra.txt
sampleFiles=$(<sra.txt)
#give 30 minutes per SRA
for sample in ${sampleFiles[*]}
	do
#	prefetch -v ${sample}
	fastq-dump --outdir . --gzip --split-files ${sample}.sra
done
rm sra.txt
