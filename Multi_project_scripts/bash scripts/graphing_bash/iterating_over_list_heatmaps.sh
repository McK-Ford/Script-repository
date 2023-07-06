#!/bin/bash
#SBATCH --mem=18G
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=deeptools_graphing_%j.out
#SBATCH --mail-type=all

module load deeptools
bd="Ordered_skew_index/"
dist_bed_list=(bydist_200.bed bydist_160.bed bydist_120.bed)
tss_bed_list=(bytss_200.bed bytss_160.bed bytss_120.bed)
pm=(plus minus)
#pm=(plus)
#classes=(200 160 120)
classes=(200)
bw_list=(skew_10 Hek_Proseq Hek_RChIP)
#bw_list=(skew_10)
for strand in ${pm[*]}
do
	for class in ${classes[*]}
		do
		for bw in ${bw_list[*]}
			do
			computeMatrix scale-regions -S bw/${bw}_${strand}.bw -R ${bd}bydist_${class}_${strand}.bed ${bd}bytss_${class}_${strand}.bed \
				-m 500 -bs 5 -b 100 -a 100 --outFileName=matrices/sorted_${bw}_${class}_${strand}.scaled.matrix

			#computeMatrix reference-point -S bw/${bw}_${strand}.bw -R ${bd}bydist_${class}_${strand}.bed ${bd}bytss_${class}_${strand}.bed \
			#-a=250 -b=250 -bs=5 --outFileName=matrices/sorted_${bw}_${class}_${strand}_TSS.matrix --referencePoint=TSS

			#computeMatrix reference-point -S bw/${bw}_${strand}.bw -R ${bd}bydist_${class}_${strand}.bed ${bd}bytss_${class}_${strand}.bed \
			#-a=250 -b=250 -bs=5 --outFileName=matrices/sorted_${bw}_${class}_${strand}_cpg_end.matrix --referencePoint=TES

			plotHeatmap -m matrices/sorted_${bw}_${class}_${strand}.scaled.matrix -out pdfs/sorted_${bw}_${class}_${strand}.scaled.pdf \
			--colorMap=Purples --sortRegions=keep
			#plotHeatmap -m matrices/sorted_${bw}_${class}_${strand}_TSS.matrix -out pdfs/sorted_${bw}_${class}_${strand}_TSS.pdf \
                        #--colorMap=Purples --sortRegions=keep --zMin 0
			#plotHeatmap -m matrices/sorted_${bw}_${class}_${strand}_cpg_end.matrix -out pdfs/sorted_${bw}_${class}_${strand}_cpg_end.pdf \
                        #--colorMap=Purples --sortRegions=keep --zMin 0
		done
	done
done

