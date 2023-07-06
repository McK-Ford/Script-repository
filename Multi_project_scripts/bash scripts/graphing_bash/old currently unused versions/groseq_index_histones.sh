#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1

module load deeptools/3.5.1

files=(G15.H3K4me2.WT G13.H3K79me2.WT G15.H3K4me3.WT)
prefix=skew_index_beds/
for file in ${files[*]}
	do
	computeMatrix scale-regions -S MCF7.Polyak.${file}.ChIPseq.dedup.rpkm.bw \
-R ${prefix}Dist_skew.bed ${prefix}TSS_skew.bed ${prefix}Low_skew_differential.bed ${prefix}neg_neutral_skew.bed -m 1000 \
-bs 10 -b 250 -a 250 --outFileName=${file}.matrix --skipZeros
	plotProfile -m ${file}.matrix -out=${file}.heatmap.pdf --outFileSortedRegions=${file}.skewindex.txt --outFileNameData=${file}.skewindex.tab.txt
done
