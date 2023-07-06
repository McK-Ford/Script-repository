#!/bin/bash

#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="import_%j.out"
#SBATCH --mail-type=all

#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588567/suppl/GSM588567_SCS814_H3K9me3_DMSO_MCF7_peak.txt.gz

#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588565/suppl/GSM588565_SCS808_H3K27me3_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588565/suppl/GSM588565_SCS808_H3K27me3_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588569/suppl/GSM588569_SCS776_H3K4me1_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588569/suppl/GSM588569_SCS776_H3K4me1_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588571/suppl/GSM588571_SCS806_H3K4me3_DMSO_806_774.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588571/suppl/GSM588571_SCS806_H3K4me3_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588573/suppl/GSM588573_SCS812_H3K9ac_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588573/suppl/GSM588573_SCS812_H3K9ac_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588575/suppl/GSM588575_SCS826_H3K14ac_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588575/suppl/GSM588575_SCS826_H3K14ac_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588575/suppl/GSM588575_SCS826_H3K14ac_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588577/suppl/GSM588577_SCS810_RNAPolII_DMSO_MCF7.bed.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM588nnn/GSM588577/suppl/GSM588577_SCS810_RNAPolII_DMSO_MCF7_peak.txt.gz
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1693nnn/GSM1693017/suppl/GSM1693017_MCF7_H3K27Ac.bw

#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM365nnn/GSM365929/suppl/GSM365929.bed.gz

#gunzip *.gz
#awk -F "\t" '{print $1, $2, $3, $5}' OFS='\t' GSM588567_SCS814_H3K9me3_DMSO_MCF7.bed > H3K9me3.bedgraph
#for f in *.bedgraph;
#	do
	#outFile=$(basename "$f" .txt | awk '{print ""$1".sorted.bedgraph"}')
	#sort -k1,1 -k2,2n $f > $outFile
#	O1=$(basename "$f" .bedgraph | awk '{print ""$1".clipped.bed"}')
#	./bedClip $f chromSizes.txt $O1
#	O2=$(basename "$f" .bedgraph | awk '{print ""$1".bw"}')
#	./bedGraphToBigWig $O1 chromSizes.txt $O2
#done
