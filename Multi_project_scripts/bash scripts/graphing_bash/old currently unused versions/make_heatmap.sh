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
regions=regions.bed
#change to relevant bed file

##############################
## Making singular heatmaps ##
##############################
#Size-sorted, edit sort order or color map if necessary. Can also add clustering when desired.
#for f in *.bw;
for sample in ${sampleFiles[*]}
	do
	#outfile1=$(basename "$f" .bw | awk '{print ""$1".matrix"}')
	echo "Generating ${sample} matrix"
	#computeMatrix reference-point -S=$f -R=CpG_genes.hg38.bed -a=3000 -b=1000 \
	#--outFileName=$outfile1 --referencePoint=TSS
	computeMatrix reference-point -S ${sample}.bw -R ${regions} -a=3000 -b=1000 \
	 --outFileName=${sample}.matrix --referencePoint=TSS
	#outfile2=$(basename "$f" .bw | awk '{print ""$1".heatmap.pdf"}')
	echo "Generating ${sample} heatmap"
	#plotHeatmap -m=$outfile1 -out=$outfile2 --colorMap=Purples --sortRegions=descend \
	# --sortUsing=region_length
	plotHeatmap -m ${sample}.matrix -out ${sample}.heatmap.pdf --colorMap=Purples \
	 --sortRegions=descend --sortUsing=region_length
done

##############################################
## Making heatmaps with shared scales/plots ##
##############################################
##If multiple heatmaps or bed files are needed

#mergeName="MCF7.H3K4me3_KDM5B" #put relevant name for merged files here

#echo "Generating ${mergeName} matrix"
#computeMatrix reference-point -S MCF7.H3K4me3.bw MCF7.KDM5B_KO.H3K4me3.bw -R=CpG_genes.hg38.bed \
# -a=3000 -b=1000 --outFileName=k4me3.matrix --referencePoint=TSS
#computeMatrix reference-point -S MCF7.H3K4me2.bw MCF7.KDM5B_KO.H3K4me2.bw -R=CpG_genes.hg38.bed \
#-a=3000 -b=1000 --outFileName=k4me2.matrix --referencePoint=TSS
#computeMatrix reference-point -S sample1.bw sample2.bw etc.bw -R region1.bed region2.bed etc.bed \
#-a=3000 -b=1000 --outFileName=${mergeName}.matrix --referencePoint=TSS

#echo "Generating ${mergeName} heatmap"
#plotHeatmap -m k4me3.matrix -out=k4me3.heatmap.pdf --colorMap=Purples --sortRegions=descend \
# --sortUsing=region_length
#plotHeatmap -m k4me2.matrix -out=k4me2.heatmap.pdf --colorMap=Purples --sortRegions=descend \
#--sortUsing=region_length
#plotHeatmap -m ${mergeName}.matrix -out=${mergeName}.heatmap.pdf --colorMap=Purples \
# --sortRegions=descend --sortUsing=region_length

###########################
## Making profile plots ##
#########################

#echo "Generating ${profile} matrix"
#computeMatrix reference-point -S sampleFiles -R CpG_genes.hg38.bed -a 3000 -b 1000 -o profile.matrix
#echo "Generating ${mergeName} heatmap"
#plotProfile -m profile.matrix --perGroup -out profile2.pdf
