#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=11162021_cpgs_proseq_noannotation_%j.out
#SBATCH --mail-type=all
module load bedtools
for i in {1..22}
do
	awk -v var="${i}" 'BEGIN{OFS='\t'; var2 = "chr"var} {if($1==var2) {print $0, $2}}' cpgIsland.hg38.txt > chr${i}_cgis.bed
	bedtools makewindows -b chr${i}_cgis.bed -w 10 -s 1 -i src > window_chr${i}.bed
	bedtools intersect -a window_chr${i}.bed -b TSS_minus.bedgraph -wa -wb > 11162021_cpgs_proseq_noanno_intersect_minus_chr${i}.txt
	bedtools intersect -a window_chr${i}.bed -b TSS_plus.bedgraph -wa -wb > 11162021_cpgs_proseq_noanno_intersect_plus_chr${i}.txt
done
#don't forget X

        awk 'BEGIN{OFS='\t'} {if($1=="chrX") print $0, $2}' cpgIsland.hg38.txt > chrX_cgis.bed
        bedtools makewindows -b chrX_cgis.bed -w 10 -s 1 -i src > window_chrX.bed
        bedtools intersect -a window_chrX.bed -b TSS_minus.bedgraph -wa -wb > 11162021_cpgs_proseq_noanno_intersect_minus_chrX.txt
        bedtools intersect -a window_chrX.bed -b TSS_plus.bedgraph -wa -wb > 11162021_cpgs_proseq_noanno_intersect_plus_chrX.txt

