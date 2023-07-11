#!/bin/bash
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=groseq_%j.out
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
module load deeptools #he uses earlier version, I will not tho
module load uni-tools/b1

# Make a sample reference sheet, lets say column 1 is file names used here
sampleFiles=$(A_Rpn5 A_Rpn9 B_Rpn6 C_Rpn7 D_Rpn8)
PE_or_SE="PE"
## Params
threads=8
umi_len=6
paired="Y"
# UMI flags, set as Y or N as appropiate
fivep_umi="Y"
threep_umi="Y"
# Adaptor sequences default to Tru-seq small RNA
adapt1 = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
adapt2= "GATCGTCGGACTGTAGAACTCTGAAC"
genome_exp = ###path to genome
genome_spike = ###path to spike in genome
#rdna = #path to rdna?
mapq=10

##############
## Pipeline ##
##############
mkdir -p logs
mkdir -p logs/fastqc

## Run fastqc on files
for file in ${sampleFiles[*]}
	do
		fastqc "${file}_R1.fastq" -o logs/fastqc
		fastqc "${file}_R2.fastq" -o logs/fastqc
	done

## Trimming adaptors and filtering rDNA reads
mkdir -p logs/fastp
mkdir -p trimmedFastq

if [[ $threep_umi == "Y" ]]
        # Branch for both UMIs
        if [[ $fivep_umi == "Y" ]]
            for pair in ${sampleFiles[*]}
                do
                    echo "trimming adapters and filtering rRNA reads for "${pair}
										fastp -i fastq/${pair}_R1.fastq -I fastq/${pair}_R2.fastq \
                        --adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
                        --umi --stdout --umi_loc=per_read --umi_len=${umi_len} \
                        --html logs/fastp/${pair}_fastp.html -w ${threads} -c \
                        --overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
                done
        # Branch for just 3' UMI
            else
            for pair in ${sampleFiles[*]}
                do
									echo "trimming adapters and filtering rRNA reads for "${pair}
									fastp -i fastq/${pair}_R1.fastq -I fastq/${pair}_R2.fastq \
										--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
										--umi --stdout --umi_loc=read1 --umi_len=${umi_len} \
										--html logs/fastp/${pair}_fastp.html -w ${threads} -c \
										--overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
                done
            fi
        # Branch for only 5' UMI or no UMIs
        else
        # Branch for only 5' UMI
        if [[ $fivep_umi == "Y" ]]
            for pair in ${sampleFiles[*]}
                do
									echo "trimming adapters and filtering rRNA reads for "${pair}
									fastp -i fastq/${pair}_R1.fastq -I fastq/${pair}_R2.fastq \
										--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
										--umi --stdout --umi_loc=read2 --umi_len=${umi_len} \
										--html logs/fastp/${pair}_fastp.html -w ${threads} -c \
										--overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
                done
                # Branch for no UMI
                else
                    for pair in ${sampleFiles[*]}
                    do
											echo "trimming adapters and filtering rRNA reads for "${pair}
											fastp -i fastq/${pair}_R1.fastq -I fastq/${pair}_R2.fastq \
												--adapter_sequence $adapt1 --adapter_sequence_r2 $adapt2 \
												--stdout \
												--html logs/fastp/${pair}_fastp.html -w ${threads} -c \
												--overlap_len_require 15 -o fastq/${pair}_R1.trimmed.fastq -O fastq/${pair}_R2.trimmed.fastq
                    done
                fi
            fi

#################################
## Mapping options, exp genome ##
#################################
mkdir -p bam
mkdir -p logs/align

if [['$PE_or_SE'=="Y"]] #is pe, Ching-hua does this step basically this same way as julius style
	for pair in ${sampleFiles[*]}
	do
			echo "aligning ${pair} to experimental genome"
			bowtie2 --local --sensitive --sensitive-local --threads ${threads} \
			-x $genome_exp -1 fastq/${pair}_R1.trimmed.fastq -2 fastq/${pair}_R2.trimmed.fastq | samtools view -bS -f 2 -q ${mapq} | samtools sort -@ ${threads} -o ${pair}.bam
			samtools index ${pair}.bam
		done
	else #PE_or_SE == "N", is SE #I like bowtie better than BWA, running that for now to keep aligner consistent. So this is more a 'julius single end' version than a 'danko single end'
		for sample in ${sampleFiles[*]}
		do
			bowtie2 --local --sensitive --sensitive-local --threads ${threads} \
			-x $genome_exp -U  fastq/${sample}_R1.trimmed.fastq | samtools view -b -q ${mapq} | samtools sort -@ ${threads} -o ${sample}.bam
			samtools index ${sample}.bam
		done
	fi

	##do spike in if I have a spikein
#	if [[ "$PAIRED" == "Y" ]]
#	    for PAIR in $(ls trimmedFastq | sed 's/_R[1-2].*//' | uniq )
#	    do
#	            echo "aligning ${PAIR} to spike in genome"
#            (bowtie2 \
#	            --local \
#	            --very-sensitive-local \
#	            --threads $(echo ${THREADS}/3*2 | bc) \
#	            --no-unal \
#	            --no-mixed \
#	            --no-discordant \
#	            -x "$GENOME_SPIKE" \
#	            -1 "trimmedFastq/${PAIR}_R1.fastq" \
#	            -2 "trimmedFastq/${PAIR}_R2.fastq" \
#	            samtools view -hS -f 2 -q ${MAPQ} |
#	            perl -n -e 'print $_ if (/^\@/ || /'${SPIKE_PREFIX}'/ ) ' |
#	            samtools view -b |
#	            samtools sort -@ $(echo ${THREADS}/3 | bc) -o spikeBAM/${PAIR}.BAM
#	            samtools index spikeBAM/${PAIR}.BAM
#	    done
#	fi
###################
## deduplication ##
###################
mkdir -p BAMdeDuped
mkdir -p logs/deDup
#mkdir -p spikeBAMdeDuped
#mkdir -p logs/spikedeDup
if [[$fivep_umi == "Y" | $threep_umi == "Y"]]
	##deduplicate with umi
	for file in ${sampleFiles[*]}
	do
		umi_tools dedup -I "$fiile" --umi-separator=":" --paired -S "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
		)> "logs/deDup/$(basename ${FILE%.BAM}_deDup.log)"
		samtools index "BAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)"
	done
else
	for file in ${sampleFiles[*]}
	do
		samtools markdup -r -@ ${threads} -o deduped.bam input.bam
		samtools index ""
	done
fi

#for FILE in spikeBAM/*.BAM
#do
#        (
#        umi_tools dedup \
#        -I "$FILE" \
#        --paired \
#        --umi-separator=":" \
#        -S "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)" \
#        )> "logs/spikedeDup/$(basename ${FILE%.BAM}_deDup.log)" &&
#        samtools index "spikeBAMdeDuped/$(basename ${FILE%.BAM}_deDuped.BAM)"
#done
#we need to skip infotable
###############################
## Making non-normalized BWs ##
###############################
mkdir -p bw/pauseloc
mkdir -p bw/10bpbins
mkdir -p bw/fulllen
for file in ${sampleFiles[*]}
	do #definitely not skipping non-covered regions, that makes it impossible to heatmap. can remove zeros later if I want.
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/pauseloc/${file}_fwd.bw --binSize 1 --numberOfProcessors ${threads} --normalizeUsing None --Offset 1 --samFlagInclude 82
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/pauseloc/${file}_rev.bw --binsize 1 --numberOfProcessors 	${threads} --normalizeUsing None --Offset 1 --samFlagInclude 98
	done
#making 10bp bins
for file in ${sampleFiles[*]}
	do
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/10bpbins/${file}_fwd_10.bw --binSize 10 --numberOfProcessors ${threads} --normalizeUsing RPKM --Offset 1 --samFlagInclude 82
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/10bpbins/${file}_rev_10.bw --binsize 10 --numberOfProcessors 	${threads} --normalizeUsing RPKM --Offset 1 --samFlagInclude 98
	done
#making full-length bigwigs
for file in ${sampleFiles[*]}
	do
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/fulllen/${file}_fwd_fl.bw --binSize 10 --numberOfProcessors ${threads} --normalizeUsing RPKM --samFlagInclude 82
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/fulllen/${file}_rev_fl.bw --binsize 10 --numberOfProcessors 	${threads} --normalizeUsing RPKM --samFlagInclude 98
	done
#rpkm norm
for file in ${sampleFiles[*]}
	do #definitely not skipping non-covered regions, that makes it impossible to heatmap. can remove zeros later if I want.
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/pauseloc/${file}_fwd.rpkm.bw --binSize 1 --numberOfProcessors ${threads} --normalizeUsing RPKM --Offset 1 --samFlagInclude 82
		bamCoverage -bam BAMdeDuped/$file.bam --outFileName bw/pauseloc/${file}_rev.rpkm.bw --binsize 1 --numberOfProcessors 	${threads} --normalizeUsing RPKM --Offset 1 --samFlagInclude 98
	done

## --ignoreDuplicates????
##Note, making the bigwigs is currently based off the julius method. The danko method is as follows:
#make bedpe using bamtools bamtobed, then use awk to select specific columns to make the beds.
#remove chromM and any rRNA
#then bedtools genomecov by strand to get bedgraphs
#then bedGraphToBigWig
