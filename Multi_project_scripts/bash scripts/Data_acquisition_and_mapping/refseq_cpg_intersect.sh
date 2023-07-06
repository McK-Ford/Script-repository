#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="hg38_cpg_gene_%j.out"
#SBATCH --mail-type=all

module load bedtools
bedtools intersect -a TSS_plus.bedgraph -b gencode_genes.hg38.txt -wa -wb > TSS_plus.txt
bedtools intersect -a TSS_minus.bedgraph -b gencode_genes.hg38.txt -wa -wb > TSS_minus.txt
#this includes all lines of all regions for every region that intersects at all with the others.
