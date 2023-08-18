#!/bin/bash
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="Acquire_Fastqs_%j.out"

#Field 2 of testref.txt is a list of SRA accession numbers. Field 1 is the desired filename. Tab separated.

module load sratoolkit
awk '{print $0=$2}' testref.txt > sra.txt
sampleFiles=$(<sra.txt)
#give 30 minutes per SRA
for sample in ${sampleFiles[*]}
	do
	prefetch -v ${sample}
	fasterq-dump --split-files ${sample} | gzip > ${sample}.fastq.gz 
done
rm sra.txt

#previous versions required use of fastq-dump or pointing to specific output directories.