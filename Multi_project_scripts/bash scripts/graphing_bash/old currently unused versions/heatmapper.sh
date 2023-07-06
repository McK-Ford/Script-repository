#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=deeptools_graphing_%j.out
#SBATCH --mail-type=all

module load deeptools/3.5.1

matrices=(G13.H3K79me2.WT G15.H3K4me2.WT G15.H3K4me3.WT skew_10 groseq_rep3 R3_rpm_proseq)
for m in ${matrices[*]}
	do
	plotHeatmap -m matrices/${m}.matrix -out ${m}.pdf --colorMap=Purples --sortRegions=descend --sortUsing=region_length
done

