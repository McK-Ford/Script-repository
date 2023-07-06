#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="skew_%j.out"
#SBATCH --mail-type=all

module load deeptools/3.5.1
computeMatrix scale-regions -S plus_strand_skew_10.bw -R TSS_CPG_end_p.hg38.bed -m 500 -bs 10 -b 250 -a 250 \
--outFileName=pos_end.matrix
computeMatrix scale-regions -S minus_strand_skew_10.bw -R TSS_CPG_end_m.hg38.bed -m 500 -bs 10 -b 250 -a 250 \
--outFileName=neg_end.matrix

#computeMatrix reference-point -S plus_strand_skew_10.bw -R TSS_CPG_end_p.hg38.bed -bs 10 -b 200 -a 4000 --outFileName=NotScaledPSkew.matrix.gz
#computeMatrix reference-point -S minus_strand_skew_10.bw -R TSS_CPG_end_m.hg38.bed -bs 10 -b 200 -a 4000 --outFileName=NotScaledMSkew.matrix.gz

computeMatrixOperations rbind -m pos_end.matrix neg_end.matrix -o end.matrix
plotProfile -m end.matrix -out=te.heatmap.pdf  --kmeans=3 \
--outFileSortedRegions=clusters.txt --silhouette --startLabel=TSS \
--endLabel=CpG_end --labelRotation=45 --samplesLabel=skew
