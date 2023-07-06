#!/bin/bash
#SBATCH --mem=256G
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output="ttseq_mapping_%j.out"
#SBATCH --mail-type=ALL

#From https://github.com/crickbabs/DRB_TT-seq, modified to work with bash/ point to the right places in my filesystem and have an all in one script.

#module load rnastar #version number important? Trying with 2.7.8a. Ching-Hua for RNA-seq uses 2.7.8a. That does appear to be the most modern version installed on bluehive. This was released feb 20 2021. Current version at master is 2.7.10b. 

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
#uncomment and modify desired method, refsheet is as usual SRRs in field 2, human IDs in field 1.
## By files in folder
#for f in *_1.fastq.gz; do echo $(basename $f _1.fastq.gz) >> "samples.txt"; done
## By refsheet
awk '{print $1}' < refSheet.txt > samples.txt

#USE FOR BOTH
sampFiles=$(<samples.txt)

echo ${sampFiles}

############
## Initialize folders and variables
###########
threads=8

ref_hg38="/scratch/mford5/references/hg38/STAR"
ref_sc3="/scratch/mford5/references/sacCer3/STAR"

rm -rf tmp

mkdir -p fastq/raw
mkdir -p fastq/trim
mkdir -p tmp
mkdir -p bam/align
mkdir -p bam/spike
mkdir -p bw
mkdir -p QC/fastp
mkdir -p QC/fragmentSize

mv *.fastq.gz fastq/raw/

######################
# get scaling factor #
######################
scaleF=1
#insert a more complicated DEseq2 r script reference here.

##################
## Make bigwigs ##
##################
echo "Making bigwigs."

for samp in ${sampFiles[*]}
do

	samtools index bam/align/${samp}.sorted.dupMark.bam

	bamCoverage --scaleFactor=${scaleF} -p ${threads} --filterRNAstrand=forward --samFlagInclude=2 --samFlagExclude=256 -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}.fwd.bw

	bamCoverage --scaleFactor=${scaleF} -p ${threads} --filterRNAstrand=forward --samFlagInclude=2 --extendReads --samFlagExclude=256 -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}_ext.fwd.bw

	bamCoverage --scaleFactor=${scaleF} -p ${threads} --filterRNAstrand=reverse --samFlagInclude=2 --samFlagExclude=256 -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}.rev.bw

	bamCoverage --scaleFactor=${scaleF} -p ${threads} --filterRNAstrand=reverse --samFlagInclude=2 --extendReads --samFlagExclude=256 -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}_ext.rev.bw
done

rm -rf tmp
