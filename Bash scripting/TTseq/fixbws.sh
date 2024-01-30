#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output="ttseq_mapping_%j.out"
#SBATCH --mail-type=ALL

#From https://github.com/crickbabs/DRB_TT-seq, modified to work with bash/ point to the right places in my filesystem and have an all in one script. Also modified to have a much more sensible normalization scheme.

#########
## Load in packages
#########

# mapping and processing
module load rnastar
module load samtools
module load picard

# qc and visualization
module load deeptools
module load fastp

# misc
module load java

###################
## Bring in data ##
###################
#uncomment and modify desired method, refsheet is as usual SRRs in field 2, human readable IDs in field 1.
## By files in folder
#for f in *_1.fastq.gz; do echo $(basename $f _1.fastq.gz) >> "samples.txt"; done
## By refsheet
awk '{print $1}' < refSheet.txt > samples.txt

#USE FOR BOTH
#declare -A sampFiles=$(<samples.txt)
mapfile -t sampFiles < samples.txt
#echo ${sampFiles}

############
## Initialize folders and variables
###########
threads=8

ref_hg38="/scratch/mford5/references/hg38/STAR"
#make a combined spike in (with chromosome prefix spike) and sample genome. Mapping separately is terrible idea as Star struggles greatly with the low coverage of the spike in.

rm -rf tmp

mkdir -p fastq/raw
mkdir -p fastq/trim
mkdir -p tmp
mkdir -p bam/align
mkdir -p bam/spike
mkdir -p bw
mkdir -p QC/fastp
mkdir -p QC/fragmentSize

mv *.fastq.gz fastq/raw

######################
# get scaling factor #
######################
module load r/4.2.1/b1
Rscript --vanilla ttseq3scale.R

mapfile -t scaleFactors < factors_tab.txt

##################
## Make bigwigs ##
##################
echo "Making bigwigs."
#mapfile -t ref2 < samples.txt
#for idx in "${!ref2[@]}"; do echo "$idx -> ${ref2[$idx]}"; done

for index in "${!sampFiles[@]}"
do
	echo "${index}"
	bm=${sampFiles[$index]}
	sf=${scaleFactors[$index]}
	echo "$bm"
	echo "$sf"
	bamCoverage --scaleFactor="$sf" -p ${threads} --filterRNAstrand=forward --samFlagInclude=2 -b bam/"${bm}".sorted.dupMark.bam -o bw/"${bm}".fwd.bw

	bamCoverage --scaleFactor="$sf" -p ${threads} --filterRNAstrand=reverse --samFlagInclude=2 -b bam/"${bm}".sorted.dupMark.bam -o bw/"${bm}".rev.bw

done
rm -rf tmp
