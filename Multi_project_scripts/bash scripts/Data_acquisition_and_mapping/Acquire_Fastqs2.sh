#!/bin/bash
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="Acquire_Fastqs_%j.out"
#SBATCH --mail-type=all

module load sratoolkit
awk '{print $0=$2}' testref.txt > sra.txt
sampleFiles=$(<sra.txt)
#give 30 minutes per SRA
for sample in ${sampleFiles[*]}
	do
	prefetch -v ${sample}
	#fastq-dump --outdir . --gzip --split-files ${sample}/${sample}.sra
	fasterq-dump --split-files ${sample} | gzip > ${sample}.fastq.gz 
done
