#!/bin/bash
#SBATCH -p standard -t 16:00:00
#SBATCH --output="CUT&Tag_mapping_%j.out"
#SBATCH -c 12 --mem=128G

#################################################################################################################
# Add --mail-type=all if want to be notified on job start and end.
# input files:
#                 1. fastq/fastq.gz
#                 2. the index reference genomes, including hg38, E. coli, and phiX
#                 3. samples.txt: contains samples' names, one per each line
#
# required programs:
#                 1. fastqc
#                 2. fastp
#                 3. bowtie2
#                 4. samtools
#                 5. picard
#                 6. deeptools
#
# note:
#                 1. bamqc needs to be installed.
#                 2. check the names of fastq/fastq.gz files.
#                 3. deeptools commands in this pipeline are based on v2.5.3. Normalized bigWig commands need
#                    to be modified if you use higher version.
#################################################################################################################

#########################################
## loading modules: required for slurm ##
#########################################

# mapping and bam file processing
module load bowtie2
module load samtools
module load picard

# QC and visualization
module load fastp
module load fastqc
module load deeptools

# misc
module load java
module load perl
module load jdk

#######################
## parameter setting ##
#######################
libraryType="PE"
numThreads=12
extendReads=200 #only used for SE #do i actually use this
bwBinSize=10
cutWindowSize=4
lengthRequired=35
mapQ=10
maxFragmentLength=800
maxin=700
# SE: single-end libs, and PE: paired-end libs; cannot use for mixed lib type

# references
reference=/scratch/mford5/references/hg38/hg38
reference_ecoli=/scratch/mford5/references/ecoli/e_coli
reference_phix=/scratch/mford5/references/phix/phix

# import the sample list
#build a list based on existing files instead ofrefsheet, uncomment
#for f in *_R1.fastq.gz; do echo $(basename $f _R1.fastq.gz) >> "samples.txt"; done
awk '{print $1}' < refSheet.txt > samples.txt
sampleFiles=$(<samples.txt)
awk '{print $0=$0".dupMark.bam"}' < samples.txt > dupMarkBamfiles.txt
awk '{print $0=$0".dedup.bam"}' < tmp_ref.txt > dedupBamfiles.txt
dupMarkBamfiles=$(<dupMarkBamfiles.txt)
dedupBamfiles=$(<dedupBamfiles.txt)

echo "Samples are:"
echo "${sampleFiles}"
echo "libraryType: ${libraryType}"
echo ""

# mkdir
mkdir -p QC/bamqc
mkdir -p QC/bowtie2_summary
mkdir -p QC/fastqc
mkdir -p QC/fastp
mkdir -p QC/fragmentSize
mkdir -p QC/misc
mkdir -p QC/multiBamSummary

mkdir -p bamFiles/dedup
mkdir -p bamFiles/dupMark
mkdir -p bamFiles/ecoli
mkdir -p bamFiles/phix

mkdir -p bigWig/dupMark/rpkm
mkdir -p bigWig/dupMark/none
mkdir -p bigWig/dedup/rpkm
mkdir -p bigWig/dedup/none

mkdir -p stats
mkdir -p rawfastqs
mkdir -p trimmedfastqs
mkdir -p tmp_files

mv *fastq.gz rawfastqs
################################
## STEP 1. QC - raw sequences ##
################################
echo "Step 1: Getting QC of raw sequences."
for FILE in rawfastqs/*fastq.gz
do
	if [ ! -s QC/fastqc/${FILE}.fastqc.html ]; then
		echo "fastqc: ${FILE}"
		fastqc "${FILE}" \
			--threads ${numThreads} \
			--outdir QC/fastqc \
			--quiet
		echo ""
	fi
done

for sample in ${sampleFiles[*]}
do
	echo "Step 1: QC and trimming ${sample} by fastp"
	if [[ "${libraryType}" == "SE" && ! -s trimmedfastqs/${sample}.trim.fastq.gz ]]
		then
		fastp -i rawfastqs/${sample}_R1.fastq.gz \
			-o trimmedfastqs/${sample}_trim.fastq.gz \
			--trim_poly_g -x --cut_window_size ${cutWindowSize} \
			--cut_tail --length_required ${lengthRequired} 2> QC/fastp/${sample}.fastp.log
#trim_poly_g removes polyG b/c in some illumina no sig rep by G
#-x removes other poly tails like polyA, cut_window_size removes low qual bases 
#in sliding windows #also check to make sure fastq names match
	elif [[ "${libraryType}" == "PE" && ! -s trimmedfastqs/${sample}_R1.trim.fastq.gz ]]
		then
		fastp --in1 rawfastqs/${sample}_R1.fastq.gz --in2 rawfastqs/${sample}_R2.fastq.gz \
			--out1 trimmedfastqs/${sample}_R1_trim.fastq.gz --out2 trimmedfastqs/${sample}_R2_trim.fastq.gz \
			--trim_poly_g -x --cut_window_size ${cutWindowSize} --cut_tail --length_required ${lengthRequired} 2> QC/fastp/${sample}.fastp.log
	else
		echo "wrong library type"
	fi 
	mv fastp.html QC/fastp/${sample}.fastp.html
	mv fastp.json QC/fastp/${sample}.fastp.json
	echo ""
done

###############################
## STEP 2. mapping - bowtie2 ##
###############################

for sample in ${sampleFiles[*]}
do
	echo "Step 2: mapping ${sample} with bowtie2"

	if [[ "${libraryType}" == "SE" && ! -s tmp_files/${sample}_fastp_trim.sam ]]
		then
		echo "mapping se to human genome"
                bowtie2 -t -p ${numThreads} -x ${reference} \
                	-U trimmedfastqs/${sample}_trim.fastq.gz -S tmp_files/${sample}_fastp_trim.sam \
                	--no-unal --local --very-sensitive-local 2> QC/bowtie2_summary/${sample}.bowtie2.txt
		echo ""
		echo "mapping to E. coli"
                bowtie2 -t -p ${numThreads} -x ${reference_ecoli} \
                	-U trimmedfastqs/${sample}_trim.fastq.gz -S tmp_files/${sample}.ecoli.sam \
                	--no-unal --local --very-sensitive-local 2> QC/bowtie2_summary/${sample}_ecoli.bowtie2.txt
		echo ""

		echo "mapping to phiX"
                bowtie2 -t -p ${numThreads} -x ${reference_phix} \
                	-U trimmedfastqs/${sample}_trim.fastq.gz -S tmp_files/${sample}.phix.sam \
                	--no-unal --local --very-sensitive-local 2> QC/bowtie2_summary/${sample}_phix.bowtie2.txt
		echo ""
	elif [[ "${libraryType}" == "PE"  && ! -s tmp_files/${sample}_fastp_trim.sam ]]
		then
		echo "mapping pe to human genome"
              bowtie2 -t -p ${numThreads} -x ${reference} \
                -1 trimmedfastqs/${sample}_R1_trim.fastq.gz -2 trimmedfastqs/${sample}_R2_trim.fastq.gz \
                -S tmp_files/${sample}.fastp_trim.sam --no-unal --local \
                --very-sensitive-local --no-mixed --no-discordant \
		--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
		echo ""
               echo "mapping to E. coli"
              bowtie2 -t -p ${numThreads} -x ${reference_ecoli} \
                        -1 trimmedfastqs/${sample}_R1_trim.fastq.gz -2 trimmedfastqs/${sample}_R2_trim.fastq.gz \
                        -S tmp_files/${sample}.ecoli.sam --no-unal --local --very-sensitive-local \
                        --no-mixed --no-discordant \
			--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
               echo ""
               echo "mapping to phiX"
               bowtie2 -t -p ${numThreads} -x ${reference_phix} \
                       -1 trimmedfastqs/${sample}_R1_trim.fastq.gz -2 trimmedfastqs/${sample}_R2_trim.fastq.gz \
                       -S tmp_files/$sample.phix.sam --no-unal --local --very-sensitive-local \
                        --no-mixed --no-discordant \
			--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
               echo ""		
	else
		echo "wrong library type or file already made"
	fi
	echo "filtering ${sample} with the mapping quality < ${mapQ}, remove unassembled contigs, mitochrondria, and chrY"
	sed '/random/d;/chrUn/d;/chrEBV/d;/chrM/d;/chrY/d' < ${sample}.fastp_trim.sam > ${sample}.filtered.sam
  	samtools view -bSq ${mapQ} tmp_files/${sample}.filtered.sam > tmp_files/${sample}.filtered.bam
	samtools view -bSq ${mapQ} tmp_files/${sample}.ecoli.sam    > tmp_files/${sample}.ecoli10.bam
	samtools view -bSq ${mapQ} tmp_files/${sample}.phix.sam     > tmp_files/${sample}.phix10.bam
	echo ""
	echo "samtools flagstat ${sample}.filtered.bam"
	samtools flagstat ${sample}.filtered.bam
	echo ""
	if [[ -s tmp_files/${sample}.filtered.bam ]]
	then
		rm -f tmp_files/${sample}.*.sam
	fi
	mv tmp_files/${sample}.ecoli10.bam bamFiles/ecoli
	mv tmp_files/${sample}.phix10.bam  bamFiles/phix
done

###########################################################################
## Step 3. bam files processing: sorting and mark duplicates with picard ##
###########################################################################

for sample in ${sampleFiles[*]}
do
	if [[ ! -s tmp_files/${sample}.sorted.bam ]]
		then
		echo "Step 3: sorting ${sample} by picard"
        	java -jar /software/picard/2.12.0/picard.jar SortSam \
        	       I=tmp_files/${sample}.filtered.bam O=tmp_files/${sample}.sorted.bam TMP_DIR=tmp \
        	       SORT_ORDER=coordinate
		rm -f tmp_files/${sample}.filtered.bam
		echo ""
	fi
	if [[ ! -s bamFiles/dupMark/${sample}.dupMark.bam ]]
	then
		echo "marking duplicates in ${sample} and indexing"
 	     	java -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
              	I=tmp_files/${sample}.sorted.bam O=bamFiles/dupMark/${sample}.dupMark.bam TMP_DIR=tmp \
              	METRICS_FILE=${sample}.dupMark.metrics use_jdk_deflater=true #see picard issue 975, script written with picard version 2.12.0
  	      	samtools index bamFiles/dupMark/${sample}.dupMark.bam
	fi
	if [[ ! -s bamFiles/dedup/${sample}.dedup.bam ]]
	then
		echo "remove duplicates in ${sample} and indexing"
		java -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
			I=tmp_files/${sample}.sorted.bam O=bamFiles/dedup/${sample}.dedup.bam TMP_DIR=tmp \
               	METRICS_FILE=${sample}.dedup.metrics REMOVE_DUPLICATES=True use_jdk_deflater=true #see picard issue 975
        	samtools index bamFiles/dedup/${sample}.dedup.bam
	fi
	echo "samtools flagstat ${sample}.dupMark.bam"
	samtools flagstat bamFiles/dupMark/${sample}.dupMark.bam
	echo ""
done

#############################################################
## Step 4: Quality controls of aligned reads (BAM files)   ##
#############################################################
echo "Step 4: multiBamSummary"
#module unload python
#module load python3
multiBamSummary bins \
	--numberOfProcessors ${numThreads} \
	--bamfiles ${dupMarkBamfiles} --labels ${sampleFiles} \
	--minmapQ ${mapQ} -o multiBamSummary.results.npz
echo ""

echo "plot heatmap for pairwise Spearman correlation coefficiencies"
plotCorrelation \
	-in multiBamSummary.results.npz --corMethod spearman \
	--skipZeros \
	--plotTitle "Spearman Correlation of Read Counts" \
	--whatToPlot heatmap --colorMap RdYlBu \
	--plotNumbers -o heatmap_SpearmanCorr_readCounts.png \
	--outFileCorMatrix SpearmanCorr_readCounts.tab
echo ""

mv multiBamSummary.results.npz         QC/multiBamSummary/
mv SpearmanCorr_readCounts.tab         QC/multiBamSummary/
mv heatmap_SpearmanCorr_readCounts.png QC/multiBamSummary/

######################################################
## Step 5: plotCoverage: coverage of the sequencing ##
######################################################
echo "Step 5: plot coverage"
echo "plotCoverage, w/ duplicates"
plotCoverage --bamfiles ${dupMarkBamfiles} --labels ${sampleFiles} --skipZeros \
        --numberOfSamples 1000000 --numberOfProcessors ${numberOfProcessors} \
        --plotFile coverage.dupMark.png --outRawCounts rawCounts.coverage.dupMark.txt
echo ""

echo "plotCoverage, w/o duplicates"
plotCoverage --bamfiles ${dedupBamfiles} --labels ${sampleFiles} --skipZeros \
        --numberOfSamples 1000000 --numberOfProcessors ${numberOfProcessors} \
        --plotFile coverage.dedup.png --outRawCounts rawCounts.coverage.dedup.txt
echo ""

mv coverage.*.png  QC/misc/
mv rawCounts.*.txt QC/misc/

##################################################################################
## Step 6: FragmentSize :the fragment size distribution of paired-end data      ##
##################################################################################
echo "Step 6: getting fragment size"
if [ "${libraryType}" == "PE" ]; then
	for sample in ${sampleFiles[*]}
	do
		bamPEFragmentSize \
			-hist fragmentSize.${sample}.png \
			--plotTitle "Fragment Size" \
			--numberOfProcessors ${numThreads} \
			--maxFragmentLength ${maxFragmentLength} \
			--bamfiles bamFiles/dupMark/${sample}.dupMark.bam \
			--samplesLabel ${sample}
		echo ""
	done
	mv fragmentSize.*.png QC/fragmentSize
fi

# plotFingerprint
plotFingerprint \
	--bamfiles bamFiles/dupMark/${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--minmapQ ${mapQ} \
	--skipZeros \
	--numberOfSamples 50000 \
	--numberOfProcessors ${numThreads} \
	--plotTitle "Fingerprints of samples" \
	--plotFile fingerprints.png \
	--outRawCounts fingerprints.tab \
	--outQualityMetrics qualityMetrics.tab

mv fingerprints.*     QC/misc/
mv qualityMetrics.tab QC/misc/

###############################################################################################################
## Step 7. coverage: to generate a sequencing-depth-normalized continuous profile of read coverages (BAM --> bigWig) ##
###############################################################################################################

echo "Step 7: making bigWig files"
echo "binSize: ${bwBinSize} and ignore chrX for normalization"
echo ""

for sample in ${sampleFiles[*]}
do
	echo "generating normalized (RPKM) bigwig files of ${sample} w/ duplicates"
	bamCoverage \
		--bam bamFiles/dupMark/${sample}.dupMark.bam \
		--outFileName bigWig/dupMark/rpkm/${sample}.dupMark.rpkm.bw \
		--binSize ${bwBinSize} \
		--extendReads \
		--numberOfProcessors ${numThreads} \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig files of ${sample} w/ duplicates"
	bamCoverage \
		--bam bamFiles/dupMark/${sample}.dupMark.bam \
		--outFileName bigWig/dupMark/none/${sample}.dupMark.bw \
		--binSize ${bwBinSize} \
		--extendReads \
		--numberOfProcessors ${numThreads}
	echo ""

	echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
	bamCoverage \
		--bam bamFiles/dedup/${sample}.dedup.bam \
		--outFileName bigWig/dedup/rpkm/${sample}.dedup.rpkm.bw \
		--ignoreDuplicates \
		--extendReads \
		--numberOfProcessors ${numThreads}
		--normalizeUsing RPKM \
		--ignoreForNormalization chrM chrX chrY
	echo ""

	echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
	bamCoverage \
		--bam bamFiles/dedup/${sample}.dedup.bam \
		--outFileName bigWig/dedup/none/${sample}.dedup.bw \
		--extendReads \
		--numberOfProcessors ${numThreads}
	echo ""
done

if [ "${libraryType}" == "SE" ]
	then
	for sample in ${sampleFiles[*]}
	do
		echo "generating normalized (RPKM) bigwig files of ${sample} w/ duplicates"
       		bamCoverage \
			--bam bamFiles/dupMark/${sample}.dupMark.bam \
			--outFileName bigWig/dupMark/rpkm/${sample}.dupMark.rpkm.bw \
                        --binSize ${bwBinSize} \
			--extendReads ${extendReads} \
                        --numberOfProcessors ${numberOfProcessors} --normalizeUsing RPKM \
                        --ignoreForNormalization chrM chrX chrY
                echo ""

		echo "generating unnormalized bigwig files of ${sample} w/ duplicates"
		bamCoverage \
			--bam bamFiles/dupMark/${sample}.dupMark.bam \
			--outFileName bigWig/dupMark/none/${sample}.dupMark.bw \
			--binSize ${bwBinSize} \
			--extendReads ${extendReads} \
			--numberOfProcessors ${numberOfProcessors}
		echo ""

		echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/dedup/${sample}.dedup.bam \
			--outFileName bigWig/dedup/rpkm/${sample}.dedup.rpkm.bw \
			--binSize ${bwBinSize} \
			--extendReads ${extendReads} \
			--numberOfProcessors ${numberOfProcessors} \
			--normalizeUsing RPKM \
			--ignoreForNormalization chrX chrM chrY
		echo ""

		echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/dedup/${sample}.dedup.bam \
			--outFileName bigWig/dedup/none/${sample}.dedup.bw \
			--binSize ${bwBinSize} \
			--extendReads ${extendReads} \
			--numberOfProcessors ${numberOfProcessors}
		echo ""
	done

elif [ "${libraryType}" == "PE" ]
	then
	for sample in ${sampleFiles[*]}
	do
		echo "generating normalized (RPKM) bigwig files of ${sample} w/ duplicates"
		bamCoverage \
			--bam bamFiles/dupMark/${sample}.dupMark.bam \
			--outFileName bigWig/dupMark/rpkm/${sample}.dupMark.rpkm.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors} \
			--normalizeUsing RPKM \
			--ignoreForNormalization chrX chrY chrM
		echo ""
		echo "generating unnormalized bigwig files of ${sample} w/ duplicates"
		bamCoverage \
			--bam bamFiles/dupMark/${sample}.dupMark.bam \
			--outFileName bigWig/dupMark/none/${sample}.dupMark.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors}
		echo ""

		echo "generating normalized (RPKM) bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/dedup/${sample}.dedup.bam \
			--outFileName bigWig/dedup/rpkm/${sample}.dedup.rpkm.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors} \
			--normalizeUsing RPKM \
			--ignoreForNormalization chrX chrY chrM
		echo ""

		echo "generating unnormalized bigwig files of ${sample} w/o duplicates"
		bamCoverage \
			--bam bamFiles/dedup/${sample}.dedup.bam \
			--outFileName bigWig/dedup/none/${sample}.dedup.bw \
			--binSize ${bwBinSize} \
			--extendReads \
			--numberOfProcessors ${numberOfProcessors}
		echo ""
	done
else
echo "wrong library type"
fi

############################################
# writing the summary for mapping pipeline #
############################################

touch stats/statsTable.tsv

echo -e \
	library'\t'\
	rawReads'\t'\
	duplicateRates'\t'\
	passedFiltersReads'\t'\
	alignConcordant'\t'\
	alignMulti'\t'\
	Unalign'\t'\
	alignConcordant%'\t'\
	alignMulti%'\t'\
	alignOverallMap%'\t' >> stats/statsTable.tsv
	
for sample in ${sampleFiles[*]}
do
	rawReads=$(cat        QC/fastp/${sample}.fastp.log              | grep "total reads:"                         | head -n 1 | awk '{print $3}')
	duplicateRate=$(cat   QC/fastp/${sample}.fastp.log              | grep "Duplication rate: "                   | head -n 1 | awk '{print $3}')
	passedReads=$(cat     QC/bowtie2_summary/${sample}.bowtie2.txt  | grep "reads; of these:$"                                | awk '{print $1}')
	alignConcordant=$(cat QC/bowtie2_summary/${sample}.bowtie2.txt  | grep "aligned concordantly exactly 1 time$"             | awk '{print $1}')
	alignMulti=$(cat      QC/bowtie2_summary/${sample}.bowtie2.txt  | grep "aligned concordantly >1 times$"                   | awk '{print $1}')
	unAlign=$(cat         QC/bowtie2_summary/${sample}.bowtie2.txt  | grep "aligned concordantly 0 times$"                    | awk '{print $1}')
	mappingRate=$(cat     QC/bowtie2_summary/${sample}.bowtie2.txt  | grep "overall alignment rate$"                          | awk '{print $1}')
	concPercent=$(echo    ${alignConcordant}"/"${passedReads}"*100" | bc -l)%
	multiPercent=$(echo   ${alignMulti}"/"${passedReads}"*100"      | bc -l)%

	echo -e \
		${sample}'\t'\
		${rawReads}'\t'\
		${duplicateRate}'\t'\
		${passedReads}'\t'\
		${alignConcordant}'\t'\
		${alignMulti}'\t'\
		${unAlign}'\t'\
		${concPercent}'\t'\
		${multiPercent}'\t'\
		${mappingRate}'\t' >> stats/statsTable.tsv
done

# multiqc to summarize fastqc, fastp, and bamqc
module unload python
module load multiqc

multiqc -d QC/fastqc -o QC -n multiqc_report_fastqc -q --no-data-dir
multiqc -d QC/fastp  -o QC -n multiqc_report_fastp  -q --no-data-dir

rm -rf tmp

echo "end of mapping"
