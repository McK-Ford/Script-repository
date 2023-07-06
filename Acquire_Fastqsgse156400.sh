#!/bin/bash
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output="Acquire_Fastqs_%j.out"
#SBATCH --mail-type=all

module load sratoolkit
sampleFiles=("SRR6904126" "SRR6904125" "SRR6904128" "SRR6904127" "SRR12467472" "SRR12467473" "SRR5364056" "SRR5364057" "SRR5364058" "SRR5364059" "SRR5364060" "SRR5364061" "SRR5364062")
#give 30 minutes per SRA
for sample in ${sampleFiles[*]}
	do
	prefetch -v ${sample}
	fasterq-dump -O . --split-files /home/mford5/ncbi/public/sra/${sample}.sra
	rm /home/mford5/ncbi/public/sra/${sample}.sra
	gzip *.fastq
done
