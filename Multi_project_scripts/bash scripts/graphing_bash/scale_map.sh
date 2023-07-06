#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1

module load deeptools/3.5.1

files=(G15.H3K4me2.WT G13.H3K79me2.WT G15.H3K4me3.WT)
prefix=skew_cluster_beds/
for file in ${files[*]}
	do
	computeMatrix scale-regions -S MCF7.Polyak.${file}.ChIPseq.dedup.rpkm.bw \
	-R ${prefix}skew_clust1.bed ${prefix}skew_clust2.bed ${prefix}skew_clust3.bed -m 1000 \
	-bs 10 -b 250 -a 250 --outFileName=${file}.matrix --skipZeros
       	plotProfile -m ${file}.matrix -out=${file}.heatmap.pdf --outFileNameData=${file}.tab.txt --startLabel=TSS \
--endLabel=CpG_end --labelRotation=45 --regionsLabel " " --plotHeight 2 --plotWidth 2

done
