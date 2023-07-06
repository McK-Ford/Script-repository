#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=deeptools_graphing_%j.out
#SBATCH --mail-type=all

module load deeptools/3.5.1
prefix=qnd_beds/qnd_
plus_minus_bws=(skew_10 groseq_rep3)
strand_indep_bws=(G15.H3K4me2.WT G13.H3K79me2.WT G15.H3K4me3.WT) 
mkdir matrices
mkdir pdfs
mkdir tabs
for file in ${plus_minus_bws[*]}
	do
	pm=(plus minus)
	#data=(bw/${file}_plus.bw bw/${file}_minus.bw)
	#regions=(${prefix}Dist_skew_p.bed ${prefix}Dist_skew_m.bed ${prefix}Dist_skew_p.bed ${prefix}Dist_skew_m.bed)
	for d in ${pm[*]}
	do
		#scaled
		computeMatrix scale-regions -S bw/${file}_${d}.bw \
-R ${prefix}Dist_skew_${d}.bed ${prefix}TSS_skew_${d}.bed -m 1000 -bs 10 -b 250 -a 250 \
--outFileName=matrices/${file}_${d}.scaled.matrix
		plotProfile -m matrices/${file}_${d}.scaled.matrix \
-out pdfs/${file}_${d}.profile.pdf --outFileNameData=tabs/${file}_${d}.scaled.tab.txt \
--startLabel=TSS --endLabel=CpG_end --labelRotation=45

#not scaled
		computeMatrix reference-point -S bw/${file}_${d}.bw -R ${prefix}Dist_skew_${d}.bed ${prefix}TSS_skew_${d}.bed -a=3000 -b=1000 \
--outFileName=matrices/${file}_${d}.matrix --referencePoint=TSS 
		plotHeatmap -m matrices/${file}_${d}.matrix -out pdfs/${file}_${d}.heatmap.pdf --colorMap=Purples \
--sortRegions=descend --sortUsing=region_length \
--outFileSortedRegions=tabs/${file}_${d}.reg.txt --outFileNameMatrix=tabs/${file}_${d}.reg.txt
	done
done

for file in ${strand_indep_bws[*]}
	do
#not scaled
	computeMatrix scale-regions -S bw/MCF7.Polyak.${file}.ChIPseq.dedup.rpkm.bw \
	-R ${prefix}Dist_skew.bed ${prefix}TSS_skew.bed -m 1000 -bs 10 -b 250 -a 250 --outFileName=matrices/${file}.matrix
	plotProfile -m matrices/${file}.matrix -out pdfs/${file}.profile.pdf \
--outFileNameData=tabs/${file}.scaled.tab.txt --startLabel=TSS --endLabel=CpG_end \
--labelRotation=45

        computeMatrix reference-point -S bw/MCF7.Polyak.${file}.ChIPseq.dedup.rpkm.bw  \
-R ${prefix}Dist_skew.bed ${prefix}TSS_skew.bed -a=3000 -b=1000 \
--outFileName=matrices/${file}.matrix --referencePoint=TSS
        plotHeatmap -m matrices/${file}.matrix -out pdfs/${file}.heatmap.pdf --colorMap=Purples \
--sortRegions=descend --sortUsing=region_length --outFileSortedRegions=tabs/${file}.reg.txt --outFileNameMatrix=tabs/${file}.reg.txt
done
