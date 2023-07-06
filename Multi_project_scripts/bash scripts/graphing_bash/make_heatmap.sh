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
sampleFiles=$(<samples.txt)

#for singular beds only
regions=bed/CpG_islands_filtered_refseq.hg38.bed

##############################
## Making singular heatmaps ##
##############################
#Size-sorted, edit sort order or color map if necessary. Can also add clustering when desired.
#for sample in *.bw;
for sample in ${sampleFiles[*]}
        do
        echo "Generating ${sample} matrix"
        computeMatrix reference-point -S ${sample} -R ${regions} -a=3000 -b=1000 \
         --outFileName=${sample}.matrix --referencePoint=TSS
        echo "Generating ${sample} heatmap"
        plotHeatmap -m ${sample}.matrix -out ${sample}.heatmap.pdf --colorMap=Purples \
         --sortRegions=descend --sortUsing=region_length
done
