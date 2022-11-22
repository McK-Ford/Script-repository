#!/bin/bash
#SBATCH -p standard
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output="proseq_mapping_QC_%j.out"
#SBATCH --mail-type=all

###############################################################################
# required files:
#                 1. fastq/fastq.gz
#                 2. the reference genome, here, hg38 is used
#                 3. annotation files: gtf and bed formats
#                 4. samples.txt: contains samples' names, one per each line
# Other concerns:
#                 1. check the names of fastq/fastq.gz files
################################################################################
module load bowtie2
module load samtools
module load deeptools
module load fastp
module load fastqc

#######################
## parameter setting ##
#######################
threads=4
spikeIn="N"
awk '{print $1}' < refSheet.txt > samples.txt
sampFiles=$(<samples.txt)
lib="SE" #PE also option

# references
ref=/scratch/mford5/references/hg38/hg38
#rdna = #path to rdna?
# Adaptor sequences default to Tru-seq small RNA
adapt1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
adapt2="GATCGTCGGACTGTAGAACTCTGAAC"

# parameters
umiLen=6
fivepUmi = "Y"
threepUmi="Y"
mapQ=10

mkdir -p fastq/trim
mkdir -p fastq/raw
mkdir -p QC/bowtie2_sum
mkdir -p QC/fastp
mkdir -p QC/fastqc
mkdir -p QC/fragmentSize
mkdir -p QC/misc
mkdir -p QC/multiBamSummary
mkdir -p bam/dedup
mkdir -p bw/none
mkdir -p bw/rpkm

mv *.fastq.gz fastq/raw
echo "Samples are:"
echo ${sampleFiles}
echo "libraryType: ${lib}"
echo ""
#################################
## STEP 1. QC1 - raw sequences ##
#################################
echo "Step 1: QC- raw sequences"
for samp in ${sampFiles[*]}
	do
	echo "trimming ${samp} by fastp and also getting fastqc"
		fastqc fastq/raw/${samp}_R1.fastq.gz -o QC/fastqc
		if [[ "${lib}" == "PE" ]]; then
			fastqc fastq/raw/${samp}_R2.fastq.gz -o QC/fastqc
		fi
	if [[ "${lib}" == "SE" && ! -s fastq/trim/${samp}_R1_trim.fastq.gz ]] ; then
		if [ ${threepUmi} == "Y" ]; then
			fastp -i fastq/raw/${samp}_R1.fastq.gz -o fastq/trim/${samp}_R1_trim.fastq.gz \
    	         		--umi --umi_len=${umiLen} --umi_loc=read1 -3 --thread ${threads} \
				--html QC/fastp/${samp}_fastp.html
		else
			fastp -i fastq/raw/${samp}_R1.fastq.gz -o fastq/trim/${samp}_R1_trim.fastq.gz \
          			-3 --thread ${threads} --html QC/fastp/${samp}_fastp.html -a ${adapt1}
		fi
	elif [["${lib}" == "PE" && ! -s fastq/trim/${samp}_R1_trim.fastq.gz ]]; then
		if [[ ${fivepUmi} == "Y" && ${threepUmi} == "Y" ]]; then
			echo "trimming adapters and filtering rRNA reads for "${samp}
			echo "using PE 5' and 3' end UMIs"
			fastp -i fastq/${samp}_R1.fastq.gz \
				-I fastq/${samp}_R2.fastq.gz \
				--adapter_sequence $adapt1 \
				--adapter_sequence_r2 $adapt2 \
				--umi --stdout --umi_loc=per_read \
				--umi_len=${umi_len} --html logs/fastp/${pair}_fastp.html \
				-w ${threads} -c --overlap_len_require 15 \
				-o fastq/${samp}_R1_trim.fastq.gz \
				-O fastq/${samp}_R2_trim.fastq.gz
		elif [[ ${fivepUmi} == "Y" && ${threepUmi} == "N" ]]; then
			echo "trimming adapters and filtering rRNA reads for "${samp}
			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${samp}_R2.fastq.gz \
				--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
				--umi --stdout --umi_loc=read2 --umi_len=${umi_len} \
				--html logs/fastp/${samp}_fastp.html -w ${threads} -c \
				--overlap_len_require 15 -o fastq/${samp}_R1_trim.fastq.gz -O fastq/${pair}_R2_trim.fastq.gz
		elif [[ ${fivepUmi} == "N" && ${threepUmi} == "Y" ]]; then
			echo "trimming adapters and filtering rRNA reads for "${samp}
			fastp -i fastq/${samp}_R1.fastq.gz -I fastq/${pair}_R2.fastq.gz \
				--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
				--umi --stdout --umi_loc=read1 --umi_len=${umi_len} \
				--html logs/fastp/${samp}_fastp.html -w ${threads} -c \
				--overlap_len_require 15 -o fastq/${samp}_R1_trim.fastq.gz -O fastq/${pair}_R2_trim.fastq.gz
		else
			echo "trimming adapters and filtering rRNA reads for "${samp}
			fastp -i fastq/${pair}_R1.fastq.gz -I fastq/${samp}_R2.fastq.gz \
				--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
				--stdout \
				--html logs/fastp/${samp}_fastp.html -w ${threads} -c \
				--overlap_len_require 15 -o fastq/${samp}_R1_trim.fastq.gz \
				-O fastq/${samp}_R2_trim.fastq.gz
	else
	echo "either wrong library type or file already trimmed"
	fi
	rm *fastp.json
    echo ""
done
###############################
## STEP 2. mapping - bowtie2 ##
###############################
#will need to edit this if I want it to handle spike in data, also chinghuas filters out rrna reads
echo "Step 2: mapping - bowtie"
for samp in ${sampFiles[*]}
do
	if [[ "${lib}" == "SE" && ! -s bam/${samp}.filtered.bam ]]; then       
		echo "mapping ${samp} to human genome"
       		bowtie2 -t -x ${ref} \
        		-U fastq/trim/${samp}_R1_trim.fastq.gz -S ${samp}.fastp_trim.sam \
        		--no-unal --local --very-sensitive-local --threads ${threads}
        	echo ""
        	echo "filtering ${samp}, remove mitochrondria and unassembled contigs"
        	sed '/chrM/d;/chrY/d;/random/d;/chrEBV/d;/alt/d;/chrUn/d' < ${samp}.fastp_trim.sam \
        	> ${samp}.filtered.sam #remove mitochondria, Y chromosome, etc.
        	samtools sort ${samp}.filtered.sam > ${samp}.sorted.sam
		samtools view -bSq ${mapQ} ${samp}.sorted.sam > bam/${samp}.filtered.bam
       		echo ""
	elif [[ "${lib}" == "PE" && ! -s bam/${samp}.filtered.bam ]]; then
		echo "mapping ${samp} to human genome"
       		bowtie2 -t -x ${ref} \
        		-1 fastq/${pair}_R1.trimmed.fastq -2 fastq/${pair}_R2.trimmed.fastq -S ${samp}.fastp_trim.sam \
        		--no-unal --local --very-sensitive-local --threads ${threads}
        	echo ""
        	echo "filtering ${samp}, remove mitochrondria and unassembled contigs"
        	sed '/chrM/d;/chrY/d;/random/d;/chrEBV/d;/alt/d;/chrUn/d' < ${samp}.fastp_trim.sam \
        	> ${samp}.filtered.sam #remove mitochondria, Y chromosome, etc.
        	samtools sort ${samp}.filtered.sam > ${samp}.sorted.sam
		samtools view -bSq ${mapQ} ${samp}.sorted.sam > bam/${samp}.filtered.bam
       		echo ""
	fi
	if [ -s bam/${samp}.filtered.bam ]
		then
       		rm -f ${samp}.*.sam
	fi
	samtools index bam/${samp}.filtered.bam
done
##################################
## STEP 3. bam files processing ##
##################################

## deduplicating with UMIs
echo "Step 3: Deduplicating with UMIs"
module unload python
module load python3
module load umi-tools/b1
for samp in ${sampFiles[*]}
do
	if [ ! -s bam/dedup/${samp}_deDuped.bam ]
	then
		umi_tools dedup \
			-I "bam/${samp}.filtered.bam" \
			--umi-separator=":" \
			-S "bam/dedup/${samp}_deDuped.bam" \
			--output-stats=${samp}_stat.txt
		samtools index bam/dedup/${samp}_deDuped.bam
	fi
done
############################################################################
## STEP 4. coverage: to generate a sequencing-depth-normalized continuous ##
## profile of read coverages (BAM --> bigWig)                             ##
############################################################################
echo "Step 7: coverage."
module unload umi-tools/b1
module load deeptools
for samp in ${sampFiles[*]}
do
	if [[ ${lib} == "SE" ]]; then
       		bamCoverage --bam bam/dedup/${samp}_deDuped.bam --outFileName bw/${samp}_fwd.bw \
                	--binSize 1 --Offset 1 --samFlagInclude 16 --numberOfProcessors ${threads} #pe 82
        	bamCoverage --bam bam/dedup/${samp}_deDuped.bam --outFileName bw/${samp}_rev.bw \
                	--binSize 1 --Offset 1 --samFlagExclude 16 --numberOfProcessors ${threads} #pe 98
	elif [[ ${lib} == "PE" ]]; then	
	        bamCoverage --bam bam/dedup/${sample}_deDuped.bam --outFileName bw/${sample}_fwd.bw \
        	        --binSize 1 --Offset 1 --samFlagInclude 82
                bamCoverage --bam bam/dedup/${sample}_deDuped.bam --outFileName bw/${sample}_rev.bw \
        	        --binSize 1 --Offset 1 --samFlagInclude 98
done
