#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output="remap_%j.out"
#SBATCH --mail-type=all

####################
## Setup + Params ##
####################
module load fastqc
module load fastp
module load bowtie2
module load samtools

awk '{print $2}' <Refsheet.txt > samples.txt
sampleFiles=$(<samples.txt)
genome="/scratch/mford5/references/hg38/hg38"

##############
## Pipeline ##
##############
mkdir -p logs
mkdir -p logs/fastqc

##Run fastqc on files
#for file in ${sampleFiles[*]}
	#do
	#	fastqc "${file}.fastq.gz" -o logs/fastqc
#done

#mkdir -p logs/fastp
#mkdir -p trimmedFastq

#for sample in ${sampleFiles[*]}
#do
	#echo "trimming adaptors in SE mode"
	#fastp -i ${sample}.fastq.gz -o ${sample}.trim.fastq.gz
#done

#############
## mapping ##
#############

mkdir -p bam
mkdir -p logs/align
#for sample in ${sampleFiles[*]}
#do 
#	echo "Aligning using SE"
#	bowtie2 --sensitive-local -x $genome -U ${sample}.trim.fastq.gz | samtools view -b | samtools sort -o ${sample}.bam
#samtools index ${sample}.bam
#done

mv *.bam bam/
mv *.bai bam/

###############################
## making non-normalized BWs ##
###############################
mkdir -p bw/sbpbins
mkdir -p bw/10bpbins

module load deeptools
plusflagInstructions=--samFlagExclude
plusflag=16
minusflag=16

for file in ${sampleFiles[*]}
do
	bamCoverage -b bam/${file}.bam -o bw/sbpbins/${file}_fwd.bw -bs 1 --normalizeUsing None --Offset 1 ${plusflagInstructions} ${plusflag}
	bamCoverage -b bam/${file}.bam -o bw/sbpbins/${file}_rev.bw -bs 1 --normalizeUsing None --Offset 1 --samFlagInclude ${minusflag}
#	bamCoverage -b bam/${file}.bam -o bw/sbpbins/${file}_fwd_10.bw -bs 10 --normalizeUsing None --Offset 1 ${plusflagInstructions} ${plusflag}
#	bamCoverage -b bam/${file}.bam -o bw/sbpbins/${file}_rev_10.bw -bs 10 --normalizeUsing None --Offset 1 --samFlagInclude ${minusflag}
done
