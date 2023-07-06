#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="too_groseq_%j.out"
#SBATCH --mail-type=all

module load deeptools/3.5.1
prefix=skew_index_beds/
computeMatrix scale-regions -S Rep3.plus.rpkm.bw -R ${prefix}Low_skew_differential_p.bed \
${prefix}neg_neutral_skew_p.bed ${prefix}Dist_skew_p.bed ${prefix}TSS_skew_p.bed -m 500 \
-bs 10 -b 250 -a 250 --outFileName=p_scaled.matrix.gz

computeMatrix scale-regions -S Rep3.minus.rpkm.bw -R ${prefix}Low_skew_differential_m.bed \
${prefix}neg_neutral_skew_m.bed ${prefix}Dist_skew_m.bed ${prefix}TSS_skew_m.bed -m 500 \
-bs 10 -b 250 -a 250 --outFileName=m_scaled.matrix.gz

plotProfile -m p_scaled.matrix.gz -out=p_scaled3.heatmap.pdf --outFileSortedRegions=groseq_p_skew.txt --outFileNameData=groseq_p_skew.tab.txt
plotProfile -m m_scaled.matrix.gz -out=m_scaled3.heatmap.pdf --outFileSortedRegions=groseq_m_skew.txt --outFileNameData=groseq_m_skew.tab.txt

