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
##############
## Initial QC
##############
#echo "QC and adapter trimming"
#for samp in ${sampFiles[*]}
#	do
#	echo "trimming ${samp} by fastp"
#	fastp -i fastq/raw/${samp}_1.fastq.gz -I fastq/raw/${samp}_2.fastq.gz \
#		-w ${threads} -c \
#		--overlap_len_require 15 -o fastq/trim/${samp}_1_trim.fastq.gz \
#		-O fastq/trim/${samp}_2_trim.fastq.gz
#	rm *fastp.json
#	mv fastp.html QC/fastp/${samp}.fastp.html
 #   echo ""
#done

##############
## Aligning ##
##############

for samp in ${sampFiles[*]}
do
	if [[ ! -s bam/align/${samp}.sorted.bam ]]
		then
		echo "mapping ${samp} to hg38"
		STAR \
			--runThreadN ${threads} --runMode alignReads --genomeDir ${ref_hg38} \
			--readFilesIn fastq/trim/${samp}_1_trim.fastq.gz fastq/trim/${samp}_2_trim.fastq.gz \
			--readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
			--twopassMode Basic --outSAMunmapped None \
			--outSAMattrRGline ID:${samp} PU:${samp} SM:${samp} LB:unknown PL:illumina \
			--outSAMtype BAM Unsorted --outTmpDir tmp/${samp} \
			--outFileNamePrefix bam/align/${samp}. \
			--sjdbGTFfile ${ref_hg38}/hg38.ensGene.gtf
		samtools sort --threads ${threads} -o bam/align/${samp}.sorted.bam bam/align/${samp}.Aligned.out.bam

		echo "mapping ${samp} to yeast"
		STAR --runThreadN ${threads} --runMode alignReads --genomeDir ${ref_sc3} \
			--readFilesIn fastq/trim/${samp}_1_trim.fastq.gz fastq/trim/${samp}_2_trim.fastq.gz \
			--readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
			--twopassMode Basic --outSAMunmapped None \
			--outSAMattrRGline ID:${samp} PU:${samp} SM:${samp} LB:unknown PL:illumina \
			--outSAMtype BAM Unsorted --outTmpDir tmp/${samp}.spike \
			--outFileNamePrefix bam/spike/${samp}. \
			--sjdbGTFfile ${ref_sc3}/sacCer3.ensGene.gtf
		samtools sort --threads ${threads} -o bam/spike/${samp}.sorted.bam bam/spike/${samp}.Aligned.out.bam

		rm bam/align/${samp}.Aligned.out.bam
		rm bam/spike/${samp}.Aligned.out.bam
	fi
done
####################
## Bam processing ##
####################

for samp in ${sampFiles[*]}
do
java -jar /software/picard/2.12.0/picard.jar MarkDuplicates INPUT=bam/align/${samp}.sorted.bam \
OUTPUT=bam/align/${samp}.sorted.dupMark.bam \
METRICS_FILE=bam/align/${samp}.sorted.marked.metrics \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=tmp/${samp}
samtools index bam/align/${samp}.sorted.marked.bam
rm bam/align/${samp}.sorted.bam

java -jar /software/picard/2.12.0/picard.jar MarkDuplicates INPUT=bam/spike/${samp}.sorted.bam \
OUTPUT=bam/spike/${samp}.sorted.dupMark.bam \
METRICS_FILE=bam/spike/${samp}.sorted.marked.metrics \
REMOVE_DUPLICATES=false ASSUME_SORTED=true \
MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=tmp/${samp}
samtools index bam/spike/${samp}.sorted.marked.bam 
rm bam/spike/${samp}.sorted.bam
done

#############################
## Fragment size distribu? ##
#############################
#Maybe don't need it but it would be nice to have to see if it compares to our other data
echo "getting fragment sizes"
for samp in ${sampFiles[*]}
do
		bamPEFragmentSize \
		-hist fragmentSize.${samp}.png \
		--plotTitle "Fragment Size" \
		--numberOfProcessors ${threads} \
		--bamfiles bam/align/${samp}.dupMark.bam \
		--samplesLabel ${samp}
	echo ""
done
mv fragmentSize.*.png QC/fragmentSize

######################
# get scaling factor #
######################
scaleF = 1
#insert a more complicated DEseq2 r script reference here.

##################
## Make bigwigs ##
##################
echo "Making bigwigs."

for samp in ${sampFiles[*]}
do

	samtools index bam/align/${samp}.sorted.dupMark.bam

	bamCoverage --scaleFactor ${scaleF} -p ${threads} --filterRNAstrand forward -b bam/align/${samp}.sorted.dupMark.bamm -o bw/${samp}.fwd.bw

	bamCoverage --scaleFactor ${scaleF} -p ${threads} --filterRNAstrand forward --extendReads -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}_ext.fwd.bw

	bamCoverage --scaleFactor ${scaleF} -p ${threads} --filterRNAstrand reverse -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}.rev.bw

	bamCoverage --scaleFactor ${scaleF} -p ${threads} --filterRNAstrand reverse --extendReads -b bam/align/${samp}.sorted.dupMark.bam -o bw/${samp}_ext.rev.bw
done

rm -rf tmp
