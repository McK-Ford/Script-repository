#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --output="skew_%j.out"
#SBATCH --mail-type=all

module load deeptools/3.5.1
computeMatrix scale-regions -S Rep3.plus.rpkm.bw \
-R skew_cluster_beds/skew_pos_clust1.bed skew_cluster_beds/skew_pos_clust2.bed skew_cluster_beds/skew_pos_clust3.bed -m 1000 \
-bs 10 -b 250 -a 250 --outFileName=pos_skew.matrix --skipZeros

computeMatrix scale-regions -S Rep3.minus.rpkm.bw \
-R skew_cluster_beds/skew_neg_clust1.bed skew_cluster_beds/skew_neg_clust2.bed skew_cluster_beds/skew_neg_clust3.bed -m 1000 \
-bs 10 -b 250 -a 250 --outFileName=neg_skew.matrix --skipZeros

computeMatrixOperations rbind -m pos_skew.matrix neg_skew.matrix -o clustered.matrix

plotProfile -m clustered.matrix -out=clustered.profile.pdf \
--outFileSortedRegions=clustered.txt --outFileNameData=clustered.tab.txt
