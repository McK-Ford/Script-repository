#!/bin/bash
#SBATCH -p standard --mem=64G  --time=8:00:00 --cpus-per-task=2
#SBATCH --output="deeptools_%j.out"

module load deeptools

peak_files=(stranded)
for file in ${peak_files[*]}
do
	computeMatrix reference-point -R gene_body_H2AZ.bed intergenic_H2AZ.bed -S ZS1_19_MCF7_EV_KDM5B.prededup.rpkm.bw ZS1_19_MCF7_EV_H2AZ.prededup.rpkm.bw ZS1_19_MCF7_EV_H3K4me3.prededup.rpkm.bw ZS1_19_MCF7_EV_H3K27me3.prededup.rpkm.bw --outFileName=${file}notscaled.jointmatrix.gz -b 2500 -a 2500 -bs 25 --referencePoint="center"
	zcat ${file}notscaled.jointmatrix.gz | sed -n "1p" > ${file}notscaledheader.txt
	zcat ${file}notscaled.jointmatrix.gz | sed "1d" | gzip > ${file}notscaled_noheader.txt.gz
	zcat ${file}notscaled_noheader.txt.gz | awk 'BEGIN{OFS="\t"}{for (i=1; i<=NF; i++) if (i<7) {printf $i"\t"} else if (i>6 && i<=206) {printf $i"\t"} else if (i>206 && i<=407) {printf $i/2"\t"} else if (i>406 && i<=607) {printf $i"\t"} else {printf $i"\t"}} {print ""}' > ${file}_modded.matrix
	cat ${file}notscaledheader.txt ${file}_modded.matrix | gzip > ${file}_notscaled_modded.matrix.gz
	plotHeatmap -m ${file}_notscaled_modded.matrix.gz --colorMap=Purples -o ${file}_nopromm_stranded_heatmap.pdf --heatmapWidth 4 --samplesLabel KDM5B H2AZ_EV H3K4me3_EV H3K27me3 --outFileSortedRegions ${file}_noprom_stranded.txt	--zMax 40
	rm ${file}notscaledheader.txt
	rm ${file}notscaled_noheader.txt.gz
	rm ${file}_modded.matrix
done
