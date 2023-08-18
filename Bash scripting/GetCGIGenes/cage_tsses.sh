#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=CpG_TSSes_%j.out

module load bedtools
for i in {1..22}
do
	awk -v var="${i}" 'BEGIN{OFS="\t"; var2 = "chr"var} {if ($1==var2) {print $1,$2,$3,$5}}' cpgIsland_hg38.txt > chr${i}_cgis.bed
	bedtools makewindows -b chr${i}_cgis.bed -w 10 -s 1 -i src > window_chr${i}.bed
	bedtools intersect -a window_chr${i}.bed -b MCF7_NET_CAGE_rev_merged.txt -wa -wb > cpgs_starts_intersect_minus_chr${i}.txt
	bedtools intersect -a window_chr${i}.bed -b MCF7_NET_CAGE_fwd_merged.txt -wa -wb > cpgs_starts_intersect_plus_chr${i}.txt
awk -v var="${i}" 'BEGIN{OFS='\t'} {if ($8!=0) {print $0}}' cpgs_starts_intersect_plus_chr${i}.txt > cpgs_starts_simple_plus_chr${i}.txt
awk -v var="${i}" 'BEGIN{OFS='\t'} {if ($8!=0) {print $0}}' cpgs_starts_intersect_minus_chr${i}.txt > cpgs_starts_simple_minus_chr${i}.txt
done

	awk 'BEGIN{OFS="\t"} {if ($1=="chrX") {print $1,$2,$3,$5}}' cpgIsland_hg38.txt > chrX_cgis.bed
        bedtools makewindows -b chrX_cgis.bed -w 10 -s 1 -i src > window_chrX.bed
        bedtools intersect -a window_chrX.bed -b MCF7_NET_CAGE_rev_merged.txt -wa -wb > cpgs_starts_intersect_minus_chrX.txt
        bedtools intersect -a window_chrX.bed -b MCF7_NET_CAGE_fwd_merged.txt -wa -wb > cpgs_starts_intersect_plus_chrX.txt

awk -v var="${i}" 'BEGIN{OFS='\t'} {if ($8!=0) {print $0}}' cpgs_starts_intersect_plus_chrX.txt > cpgs_starts_simple_plus_chrX.txt
awk -v var="${i}" 'BEGIN{OFS='\t'} {if ($8!=0) {print $0}}' cpgs_starts_intersect_minus_chrX.txt > cpgs_starts_simple_minus_chrX.txt

rm *intersect*
