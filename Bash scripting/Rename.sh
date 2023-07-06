#!/bin/bash
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="rename_%j.out"
#SBATCH --mail-type=all

#This file is useful for taking a tab-separated 
#reference sheet with filename prefixes (ex SRR version & human readable version)
#and converting the related files over.

while read outFile inFile; do mv "${inFile}_1.fastq.gz" "${outFile}.fastq.gz"; done < Ref_sheet_2.txt

