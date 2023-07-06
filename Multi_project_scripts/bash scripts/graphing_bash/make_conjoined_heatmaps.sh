#!/bin/bash

#SBATCH --mem=12G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="heatmaps_%j.out"
#SBATCH --mail-type=all

###############################
## Loading necessary modules ##
###############################
module load deeptools

#####################
## Getting samples ##
#####################
#for singular heatmaps / profile map only
#sampleFiles=$(<samples.txt)
#sampleFiles= ls *.bw
#for singular beds only
#regions=TSS_CpG_end_refseq.hg38.bed
#change to relevant bed file

##############################################
## Making heatmaps with shared scales/plots ##
##############################################
##If multiple heatmaps or bed files are needed

#mergeName="MCF7.H3K4me3_KDM5B" #put relevant name for merged files here

#echo "Generating ${mergeName} matrix"
computeMatrix reference-point \
-S MCF7.Polyak.G15.H3K4me3.WT.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N1.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me3.KDM5B_KO_N2.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_0d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_2d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me3.C70_7d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3k4me3.C70_14d.ChIPseq.dedup.rpkm.bw \
-R=../bed/TSS_CpG_end_refseq.hg38.bed \
-a=3000 -b=1000 --outFileName=h3k4me3.matrix --referencePoint=TSS

computeMatrix reference-point \
-S MCF7.Polyak.G15.H3K4me2.WT.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me2.KDM5B_KO_N1.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G15.H3K4me2.KDM5B_KO_N2.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_0d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_2d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_7d.ChIPseq.dedup.rpkm.bw \
MCF7.Polyak.G1.H3K4me2.C70_14d.ChIPseq.dedup.rpkm.bw \
-R=../bed/TSS_CpG_end_refseq.hg38.bed \
-a=3000 -b=1000 --outFileName=h3k4me2.matrix --referencePoint=TSS

plotHeatmap -m h3k4me3.matrix -out=h3k4me3.heatmap.pdf --colorMap=Purples --sortRegions=descend \
--sortUsing=region_length --samplesLabel H3K4me3_WT H3K4me3_KDM5B_KO1 H3K4me3_KDM5B_KO1 H3K4me3_C70_0 H3K4me3_C70_2 H3K4me3_C70_7 H3K4me3_C70_14

plotHeatmap -m h3k4me2.matrix -out=h3k4me2.heatmap.pdf --colorMap=Purples --sortRegions=descend \
--sortUsing=region_length --samplesLabel H3K4me2_WT H3K4me2_KDM5B_KO1 H3K4me2_KDM5B_KO2 H3K4me2_C70_0 H3K4me2_C70_2 H3K4me2_C70_7 H3K4me2_C70_14

