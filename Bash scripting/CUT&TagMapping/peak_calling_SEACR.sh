#!/bin/bash
#SBATCH -p standard --mem=32G  --time=1:00:00 --cpus-per-task=1
#SBATCH --output="SEACR_peak_calling_%j.out"

module load seacr
module load R

files=(EV_KDM5B)
for file in ${files[*]}
do
	tools/bigWigToBedGraph ZS1_19_MCF7_${file}.prededup.rpkm.bw ${file}.bedGraph
	awk '{if ($4!=0) {print $0}}' ${file}.bedGraph > ${file}.no0.bedGraph
	SEACR_1.3.sh ${file}.no0.bedGraph 0.0025 non stringent ZS1_19_MCF7_0025${file}
	rm *.bedGraph
done
