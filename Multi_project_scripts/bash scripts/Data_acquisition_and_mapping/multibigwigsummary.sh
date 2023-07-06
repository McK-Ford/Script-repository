#!/bin/bash
#SBATCH -p standard
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output="multibigwigsummary.out"
#SBATCH --mail-type=all

multiBigwigSummary bins -b ZS1_19_MCF7_EV_KDM5B.dedup.rpkm.bw ZS1_19_MCF7_EV_H2AZ.dedup.rpkm.bw ZS1_19_MCF7_EV_H3K4me3.dedup.rpkm.bw ZS1_19_MCF7_SH1_H2AZ.dedup.rpkm.bw ZS1_19_MCF7_SH1_H3K4me3.dedup.rpkm.bw ZS1_19_MCF7_SH2_H2AZ.dedup.rpkm.bw ZS1_19_MCF7_SH2_H3K4me3.dedup.rpkm.bw -l EV_KDM5B EV_H2AZ EV_H3K4me3 S1_H2AZ SH1_H3K4me3 SH2_H2AZ SH2_H3K4me3 -o rpkmbws.npz -p max

