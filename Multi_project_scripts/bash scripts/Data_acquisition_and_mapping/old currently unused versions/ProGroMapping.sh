#!/bin/bash
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=proseq_%j.out
#SBATCH --mail-type=all

####################
## Setup + Params ##
####################
## Loading in modules
## For Ching-hua heavily julius based version
module load fastqc
module load fastp
module load bowtie2
module load samtools
module load umi-tools/b1

# Make a sample reference sheet, lets say column 1 is file names used here
awk '{print $1}' < DankoRefSheet.txt > samples.txt
sampleFiles=$(<samples.txt)
PE_or_SE="SE"
## Params
threads=8
umi_len=6
paired="N"
# UMI flags, set as Y or N as appropiate
fivep_umi="N"
threep_umi="N"
# Adaptor sequences default to Tru-seq small RNA
adapt1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
adapt2="GATCGTCGGACTGTAGAACTCTGAAC"
genome_exp="/scratch/mford5/references/hg38/hg38"
mapq=10

##############
## Pipeline ##
##############
mkdir -p logs
mkdir -p logs/fastqc

## Run fastqc on files
#for file in ${sampleFiles[*]}
#	do
#		fastqc "${file}_R1.fastq.gz" -o logs/fastqc
#	if [[ ${PE_or_SE} == "PE" ]]
#	then
#		fastqc "${file}_R2.fastq.gz" -o logs/fastqc
#	fi
#done

## Trimming adaptors and filtering rDNA reads
mkdir -p logs/fastp
mkdir -p trimmedFastq

#if [[ $threep_umi == "Y" ]]
        # Branch for both UMIs
#        then
#	if [[ $fivep_umi == "Y" ]]
#		then
#		for pair in ${sampleFiles[*]}
 #           do
#			echo "trimming adapters and filtering rRNA reads for "${pair}
#			echo "using PE 5' and 3' end UMIs"
#			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${pair}_R2.fastq.gz --adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 --umi --stdout --umi_loc=per_read --umi_len=${umi_len} --html logs/fastp/${pair}_fastp.html -w ${threads} -c --overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
#                done
        # Branch for just 3'
#	else
#		for pair in ${sampleFiles[*]}
 #               do
	#		echo "trimming adapters and filtering rRNA reads for "${pair}
#			echo "using PE 3' only UMI"	
#			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${pair}_R2.fastq.gz --adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 --umi --stdout --umi_loc=read1 --umi_len=${umi_len} --html logs/fastp/${pair}_fastp.html -w ${threads} -c --overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
#                done
#        fi
        # Branch for only 5' UMI or no UMIs
#else
        # Branch for only 5' UMI
#        if [[ $fivep_umi == "Y" ]]
#		then
 #           for pair in ${sampleFiles[*]}
 #               do
#			echo "trimming adapters and filtering rRNA reads for "${pair}
#			echo "Using PE 5' only UMI"
#			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${pair}_R2.fastq.gz --adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 --umi --stdout --umi_loc=read2 --umi_len=${umi_len} --html logs/fastp/${pair}_fastp.html -w ${threads} -c --overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
#                done
                # Branch for no UMI
#           else
#		if [[ ${PE_or_SE} == "PE" ]]
#		then
#			for pair in ${sampleFiles[*]}
#			do
#			echo "trimming adapters and filtering rRNA reads for "${pair}
#			echo "Using PE no UMI"
#			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${pair}_R2.fastq.gz --adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 --html logs/fastp/${pair}_fastp.html -w ${threads} -c --overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
#			done
#		else
#			for sample in ${sampleFiles[*]}
#			do
#			echo "Using SE"
#			fastp -i ${sample}_R1.fastq.gz -o ${sample}.trim.fastq.gz
#			done
#		fi                    
 #         fi
#fi


#################################
## Mapping options, exp genome ##
#################################
mkdir -p bam
mkdir -p logs/align

#if [[ $PE_or_SE == "PE" ]] #is pe, Ching-hua does this step basically same way as julius style
#	then
#	for pair in ${sampleFiles[*]}
#	do
#		echo "aligning ${pair} to experimental genome using PE"
#		bowtie2 --sensitive-local --threads ${threads} -x $genome_exp -1 fastq/${pair}_R1.trimmed.fastq -2 fastq/${pair}_R2.trimmed.fastq | samtools view -bS -f 2 -q ${mapq} | samtools sort -@ ${threads} -o ${pair}.bam
#		samtools index ${pair}.bam
#	done
#else #PE_or_SE == "N", is SE #I like bowtie better than BWA, want to keep aligner consistent. So this is more a 'julius single end' version than a 'danko single end'
#	for sample in ${sampleFiles[*]}
#	do
#		echo "Aligning using SE"
#		bowtie2 --sensitive-local --threads ${threads} -x $genome_exp -U ${sample}_R1.trim.fastq.gz | samtools view -b -q ${mapq} | samtools sort -@ ${threads} -o ${sample}.bam
#		samtools index ${sample}.bam
#	done
#
mv *.bam bam/
mv *.bai bam/
###################
## deduplication ##
###################
mkdir -p bamdeduped
mkdir -p logs/deDup
#if [[$fivep_umi == "Y" | $threep_umi == "Y"]]
#	then
	##deduplicate with umi
#	for file in ${sampleFiles[*]}
#	do
#		umi_tools dedup -I "$file" --umi-separator=":" --paired -S "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)"
#		samtools index "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)"
#	done
#else
#	for file in ${sampleFiles[*]}
#	do	
#		samtools sort -n -o bam/${file}.namesort.bam bam/${file}.bam
#		samtools fixmate -m bam/${file}.namesort.bam bam/${file}.wmates.bam
#		samtools sort -o bam/${file}.possort.bam bam/${file}.wmates.bam
#		rm bam/${file}.namesort.bam bam/${file}.wmates.bam
#		samtools markdup bam/${file}.bam bamdeduped/${file}.deduped.bam -r -@ ${threads}
#		samtools index bamdeduped/${file}.deduped.bam
#	done
#fi

###############################
## Making non-normalized BWs ##
###############################
mkdir -p bw/pauseloc
mkdir -p bw/10bpbins

module unload umi-tools
module load deeptools
if [[ ${PE_or_SE} == "PE" ]]
	then
	plusflagInstructions=--samFlagInclude
	plusflag=82
	minusflag=98
else
	plusflagInstructions=--samFlagExclude
	plusflag=16
	minusflag=16
fi

for file in ${sampleFiles[*]}
	do #definitely not skipping non-covered regions, that makes it impossible to heatmap. can remove zeros later if I want.
		bamCoverage -b bam/${file}.bam -o bw/pauseloc/${file}_fwd.bw -bs 1 --numberOfProcessors ${threads} --normalizeUsing None --Offset 1 ${plusflagInstructions} ${plusflag}
		bamCoverage -b bam/${file}.bam -o bw/pauseloc/${file}_rev.bw -bs 1 --numberOfProcessors 	${threads} --normalizeUsing None --Offset 1 --samFlagInclude ${minusflag}
#making 10bp bins
		bamCoverage -b bam/${file}.bam -o bw/10bpbins/${file}_fwd_10.bw -bs 10 --numberOfProcessors ${threads} --normalizeUsing None --Offset 1 ${plusflagInstructions} ${plusflag}
		bamCoverage -b bam/${file}.bam -o bw/10bpbins/${file}_rev_10.bw -bs 10 --numberOfProcessors 	${threads} --normalizeUsing None --Offset 1 --samFlagInclude ${minusflag}
	done
