#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output="Acquire_Fastqs_%j.out"


module load sratoolkit
awk '{print $2}' < refSheet.txt > tmp_SRR.txt
sampleFiles=$(<tmp_SRR.txt)
for sample in ${sampleFiles[*]}
	do
	fasterq-dump --split-files ${sample}
done
