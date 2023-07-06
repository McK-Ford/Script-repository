#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=deeptools_graphing_%j.out
#SBATCH --mail-type=all

module load deeptools
bed_prefix="08192021_01bidirectional_islands/"
bw_prefix="bw/groseq_rep3"
plus_minus_beds=(bidir unidir close_bidir far_bidir)
pm=(plus minus)
for d in ${pm[*]}
do
#scaled
	#sort-bed ${bed_prefix}bidir_${d}.bed >${bed_prefix}bidir_${d}.sorted.bed
	#sort-bed ${bed_prefix}unidir_${d}.bed >${bed_prefix}unidir_${d}.sorted.bed
	#sort-bed ${bed_prefix}close_bidir_${d}.bed >${bed_prefix}close_bidir_${d}.sorted.bed
	#sort-bed ${bed_prefix}far_bidir_${d}.bed >${bed_prefix}far_bidir_${d}.sorted.bed
#	computeMatrix scale-regions -S ${bw_prefix}_${d}.bw \
#	-R ${bed_prefix}bidir_${d}.sorted.bed ${bed_prefix}unidir_${d}.sorted.bed ${bed_prefix}close_bidir_${d}.sorted.bed ${bed_prefix}far_bidir_${d}.sorted.bed \
#	-m 1000 -bs 10 -b 250 -a 250 \
#	--outFileName=matrices/bidir_islands_groseq_${d}.scaled.matrix
	
#	plotProfile -m matrices/bidir_islands_groseq_${d}.scaled.matrix \
#	-out pdfs/bidir_islands_groseq_${d}.profile.pdf --outFileNameData=tabs/bidir_islands_groseq_${d}.scaled.tab.txt \
#	--startLabel=TSS --endLabel=CpG_end --labelRotation=45

#not scaled
#	computeMatrix reference-point -S ${bw_prefix}_${d}.bw -R ${bed_prefix}bidir_${d}.sorted.bed ${bed_prefix}unidir_${d}.sorted.bed \
#	${bed_prefix}close_bidir_${d}.sorted.bed ${bed_prefix}far_bidir_${d}.sorted.bed -a=3000 -b=1000 \
#	--outFileName=matrices/bidir_islands_groseq_${d}.matrix --referencePoint=TSS 
	
	plotHeatmap -m matrices/bidir_islands_groseq_${d}.matrix -out pdfs/bidir_islands_groseq_${d}.heatmap.pdf --colorMap=Purples \
	--sortRegions=descend --sortUsing=region_length \
	--outFileSortedRegions=tabs/bidir_islands_groseq_${d}.reg.txt --outFileNameMatrix=tabs/bidir_islands_groseq_${d}.reg.txt --zMax 10
done

