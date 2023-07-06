#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="skew_%j.out"
#SBATCH --mail-type=all

module load deeptools
#computeMatrix reference-point -S pos_skew.bw -R plus_only.bed -a=3000 -b=1000 --outFileName=plus_skew.matrix --referencePoint=TSS
#plotHeatmap -m plus_skew.matrix -out=plus_skew.heatmap.pdf --colorMap=Purples --sortRegions=descend --sortUsing=region_length --kmeans=3

#computeMatrix reference-point -S neg_skew.bw -R minus_only.bed -a=3000 -b=1000 --outFileName=minus_skew.matrix --referencePoint=TSS
#plotHeatmap -m minus_skew.matrix -out=minus_skew.heatmap.pdf --colorMap=Purples --sortRegions=descend --sortUsing=region_length --kmeans=3

computeMatrix reference-point -S senseskew_small.bw -R bed/TSS_CpG_end_refseq.hg38.bed -a=1000 -b=250 --outFileName=senseskew_small.matrix --referencePoint=TSS
plotHeatmap -m senseskew_small.matrix -out=senseskew_small.heatmap.pdf --colorMap=Purples --sortRegions=descend --sortUsing=region_length --kmeans=3

computeMatrix reference-point -S senseskew_small.bw -R bed/TSS_CpG_end_refseq.hg38.bed -a=3000 -b=1000 --outFileName=senseskew_small1.matrix --referencePoint=TSS
plotHeatmap -m senseskew1_small.matrix -out=senseskew_small1.heatmap.pdf --colorMap=Purples --sortRegions=descend --sortUsing=region_length

#computeMatrix scale-regions -S senseskew.bw -R bed/TSS_CpG_end_refseq.hg38.bed --outFileName=scaled_sense.matrix -b=250 -a=250 -m=750
#plotHeatmap -m scaled_sense.matrix -out=sense_scaled_clustered.heatmap.pdf --colorMap=Purples  --kmeans=3

