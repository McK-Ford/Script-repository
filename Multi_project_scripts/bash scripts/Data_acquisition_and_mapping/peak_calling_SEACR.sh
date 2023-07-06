#!/bin/bash
#SBATCH -p standard
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output="SEACR_peak_calling_%j.out"
#SBATCH --mail-type=all

module load seacr
module load R

files=(EV_H2AZ SH1_H2AZ SH2_H2AZ EV_H3K4me3 SH1_H3K4me3 SH2_H3K4me3)
for file in ${files[*]}
do
	../tools/bigWigToBedGraph ZS1_19_MCF7_${file}.dedup.rpkm.bw ZS1_19_MCF7_${file}.bedGraph #seacr requires bedgraph
	awk '{if ($4!=0) {print $0}}' ZS1_19_MCF7_${file}.bedGraph > ZS1_19_MCF7_${file}.NO_0.bedGraph #seacr needs the bedgraph to not have 0 scores, just leave those spots empty
	SEACR_1.3.sh ZS1_19_MCF7_${file}.NO_0.bedGraph 0.01 non stringent ZS1_19_MCF7_${file}
done
