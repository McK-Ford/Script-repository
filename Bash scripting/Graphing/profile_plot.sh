#!/bin/bash
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="profile_%j.out"
#SBATCH --mail-type=all

module load deeptools

computeMatrix scale-regions \
-S MCF7.Polyak.G15.H3K4me3.WT.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N1.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N2.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_0d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3k4me3.C70_14d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_2d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_7d.ChIPseq.dedup.rpkm.bw \
-R=../bed/TSS_CpG_end_refseq.hg38.bed \
-m 500 -bs 10 -b 250 -a 250 -o profile3.matrix

computeMatrix scale-regions \
-S MCF7.Polyak.G15.H3K4me2.WT.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me2.KDM5B_KO_N1.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me2.KDM5B_KO_N2.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_0d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_14d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_2d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_7d.ChIPseq.dedup.rpkm.bw \
-R=../bed/TSS_CpG_end_refseq.hg38.bed \
-m 500 -bs 10 -b 250 -a 250 -o profile2.matrix


plotProfile -m profile2.matrix --perGroup -out profile2.pdf \
--samplesLabel H3K4me2_WT H3K4me2_N1_KDM5B_KO H3K4me2_N2_KDM5B_KO H3K4me2_C70_0 \
H3K4me2_C70_14 H3K4me2_C70_2 H3K4me2_C70_7 --startLabel=TSS --endLabel=CpG_end

plotProfile -m profile3.matrix --perGroup -out profile3.pdf \
--samplesLabel H3K4me3_WT H3K4me3_N1_KDM5B_KO H3K4me3_N2_KDM5B_KO H3K4me3_C70_0 \
H3K4me3_C70_14 H3K4me3_C70_2 H3K4me3_C70_7 --startLabel=TSS --endLabel=CpG_end
