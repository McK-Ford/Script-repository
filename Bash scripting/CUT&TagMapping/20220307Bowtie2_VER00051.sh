#!/bin/bash
#SBATCH -p standard
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output="Bowtie_mapping_QC_%j.out"
#SBATCH --mail-type=all

###################################################################################################
# required files:
#                 1. fastq/fastq.gz
#                 2. the reference genome, here, hg38 is used
#                 3. annotation files: gtf and bed formats
#                 4. samples.txt: contains samples' names, one per each line
# Notice:
#                 1. fastp needs to be installed - check to make sure its location relative to
#                      script is correct.
#                 2. check the names of fastq/fastq.gz files; the file names in this scripts is the
#                    output from URMC genome sequencing core.
####################################################################################################
#multiply time by number of jobs.
module load bowtie2
module load samtools
module load picard
module load deeptools
module load java
module load fastp
module load perl

#######################
## parameter setting ##
#######################
# library
libraryType="PE"
# SE: single-end libraries, and PE: paired-end libraries; cannot be used for mixed libraries
spikeIn="NO"
extendReads=200 #only used in SE libraries.
awk '{print $2}' > tmp_ref.txt < sampleRefSheet.txt
sampleFiles=$(<tmp_ref.txt)

# references
reference=/scratch/mford5/references/hg38/hg38
reference_ecoli=/scratch/mford5/references/ecoli/e_coli
reference_phix=/scratch/mford5/references/phix/phix

# parameters
bwBinSize=10
mappingQuality=10
maxin=700
numberOfProcessors=16
lengthRequired=35
maxFragmentLength=800
cutWindowSize=4

# import the sample list
awk '{print $0=$0".bam"}' < sampleRefSheet.txt > bamfiles.txt #Hmm... do I actually need this
awk '{print $0=$0".dupMark.bam"}' < sampleRefSheet.txt > dupMarkBamfiles.txt
awk '{print $0=$0".dedup.bam"}' < sampleRefSheet.txt > dedupBamfiles.txt


dupMarkBamfiles=$(<dupMarkBamfiles.txt)
dedupBamfiles=$(<dedupBamfiles.txt)


mkdir fastq
mkdir -p QC/bowtie2_summary
mkdir -p QC/fastp
mkdir -p QC/fragmentSize
mkdir -p QC/misc
mkdir -p QC/multiBamSummary
mkdir -p bamFiles/dedup
mkdir -p bamFiles/dupMark
mkdir -p bamFiles/ecoli
mkdir -p bamFiles/phix
mkdir -p bigWig/dupMark/rpkam
mkdir -p bigWig/dupMark/none
mkdir -p bigWig/dedup/rpkm
mkdir -p bigWig/dedup/none
mkdir stats
mkdir raw
mv *.fastq.gz fastq

echo "Samples are:"
echo ${sampleFiles}
echo "libraryType: ${libraryType}"
echo "bwBinSize = ${bwBinSize}"
echo "maximun fragment length: ${maxin}"
echo "mapping quality cut-off: ${mappingQuality}"
echo ""
echo ""
#################################
## STEP 1. QC1 - raw sequences ##
#################################
#with 16 gb and 4 CPU this step takes ~ 1 hr for 13 (miseq) samples.
echo "Step 1: QC- raw sequences"
for sample in ${sampleFiles[*]}
	do
	echo "trimming ${sample} by fastp"
	if [ "${libraryType}" == "SE" ]; then
		/scratch/mford5/tools/fastp \
			-i fastq/${sample}.fastq.gz -o fastq/${sample}.trim.fastq.gz \
			--trim_poly_g -x --cut_window_size ${cutWindowSize} --cut_tail --length_required 35

        elif [ "${libraryType}" == "PE" ]; then
		/scratch/mford5/tools/fastp --in1 fastq/${sample}_R1.fastq.gz \
			--in2 fastq/${sample}_R2.fastq.gz --out1 fastq/${sample}_R1_trim.fastq.gz \
                        --out2 fastq/${sample}_R2_trim.fastq.gz --trim_poly_g -x \
			--cut_window_size ${cutWindowSize} \
                        --cut_tail --length_required ${lengthRequired} 2> QC/fastp/${sample}.fastp.log
        else
               echo "wrong library type"
        fi

	mv fastp.html ${sample}.fastp.html
        mv fastp.json ${sample}.fastp.json
        echo ""
done
mv *.fastp.html QC/fastp
mv*.fastp.json QC/fastp
##Murphy lab does not include this trimming step, because bowtie trims many standard adapters on its own.
###############################
## STEP 2. mapping - bowtie2 ##
###############################

echo "Step 2: mapping - bowtie"
for sample in ${sampleFiles[*]}
do
	echo "mapping ${sample}" #announces name of sample

##THIS SECTION HANDLES SE LIBRARIES
        if [ "${libraryType}" == "SE" ]; then          
		echo "mapping ${sample} to human genome"
                bowtie2 -t -p ${numberOfProcessors} -x ${reference} \
                	-U fastq/${sample}.trim.fastq.gz -S ${sample}.fastp_trim.sam \
                	--no-unal --local --very-sensitive-local
               	echo ""
                echo "mapping ${sample} to E. coli"
                bowtie2 -t -p ${numberOfProcessors} -x ${reference_ecoli} \
                	-U fastq/${sample}.trim.fastq.gz      -S ${sample}.ecoli.sam \
                	--no-unal --local --very-sensitive-local
                echo ""

                echo "mapping ${sample} to phix"
                bowtie2 -t -p ${numberOfProcessors} -x ${reference_phix} \
                	-U fastq/${sample}.trim.fastq.gz -S ${sample}.phix.sam \
                	--no-unal --local --very-sensitive-local
                echo ""

#-x reference genome, -U is input, -S is output, -t is time. -no-unal is no unaligned written to
# the output sam file. --local just makes it so it's local alignment as expected. This means reads
# can be trimmed for alignment instead of needing exact end matches.

        elif [ "${libraryType}" == "PE" ]; then
               echo "mapping ${sample} to hg"
              bowtie2 -t -p ${numberOfProcessors} -x ${reference} \
                -1 fastq/${sample}_R1_trim.fastq.gz -2 fastq/${sample}_R2_trim.fastq.gz \
                -S ${sample}.fastp_trim.sam --no-unal --local \
                --very-sensitive-local --no-mixed --no-discordant \
		--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
            echo ""
#so we add the sam output, the no unal/mixed/discordant arguments, the very sensitive local argument, and the limitation on
#maximum fragment length to patrick's code. It's unclear however if he does any of this when he pipes to SAM because that's
#not in the screenshot.

               echo "mapping ${sample} to E. coli"
              bowtie2 -t -p ${numberOfProcessors} -x ${reference_ecoli} \
                        -1 fastq/${sample}_R1_trim.fastq.gz -2 fastq/${sample}_R2_trim.fastq.gz \
                       -S ${sample}.ecoli.sam --no-unal --local --very-sensitive-local \
                        --no-mixed --no-discordant \
			--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
                echo ""
                echo "mapping ${sample} to phiX"
                bowtie2 -t -p ${numberOfProcessors} -x ${reference_phix} \
                        -1 fastq/${sample}_R1_trim.fastq.gz -2 fastq/${sample}_R2_trim.fastq.gz \
                       -S $sample.phix.sam --no-unal --local --very-sensitive-local \
                        --no-mixed --no-discordant \
			--maxins ${maxin} 2> QC/bowtie2_summary/${sample}.bowtie2.txt
                echo ""
        else
                echo "wrong library type"
       fi

##maxins is the maximum fragment length. -no-discordant is nothing aligns singly, and no mixed too.

        echo "filtering ${sample} with mapping quality < ${mappingQuality}, remove mitochrondria and unassembled contigs"
        sed '/chrM/d;/chrY/d;/random/d;/chrEBV/d;/alt/d;/chrUn/d' < ${sample}.fastp_trim.sam \
        > ${sample}.filtered.sam #remove mitochondria, Y chromosome, etc.
        samtools view -bSq ${mappingQuality} ${sample}.filtered.sam > ${sample}.filtered.bam
        samtools view -bSq ${mappingQuality} ${sample}.ecoli.sam    > ${sample}.ecoli10.bam
        samtools view -bSq ${mappingQuality} ${sample}.phix.sam     > ${sample}.phix10.bam
        echo ""
# -bSq is output in bam format(b), skip alignments where mapq is low (q), and S is compatibility.
        rm -f ${sample}.*.sam
done
##################################
## STEP 3. bam files processing ##
##################################

echo "Step 3: bam file processing"
for sample in ${sampleFiles[*]}
	do
        echo "sorting ${sample} by picard"
        java -jar /software/picard/2.12.0/picard.jar SortSam \
               I=${sample}.filtered.bam O=${sample}.sorted.bam TMP_DIR=tmp \
               SORT_ORDER=coordinate
        rm -f ${sample}.filtered.bam
        echo ""
#
        echo "marking duplicates in ${sample} and indexing"
        java -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
               I=${sample}.sorted.bam O=${sample}.dupMark.bam TMP_DIR=tmp \
               METRICS_FILE=${sample}.dupMark.metrics
        samtools index ${sample}.dupMark.bam
        echo ""

	echo "remove duplicates in ${sample} and indexing" #see deeptools issue 955 for why remove duplicates now instead of when making bigwig.
	java -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
		I=${sample}.sorted.bam O=${sample}.dedup.bam TMP_DIR=tmp \
               METRICS_FILE=${sample}.dedup.metrics REMOVE_DUPLICATES=True
        samtools index ${sample}.dedup.bam
        rm -f ${sample}.sorted.bam
        echo ""

	samtools flagstat ${sample}.dupMark.bam
done

##############################################################################################
## STEP 4. multiBamSummary: similarity of read distribution in different sequencing samples ##
##############################################################################################
echo "Step 4: multiBamSummary and plot correlation heatmap"
multiBamSummary bins --numberOfProcessors ${numberOfProcessors} \
        --bamfiles ${dupMarkBamfiles} --labels ${sampleFiles} \
        --minMappingQuality ${mappingQuality} -o multiBamSummary.results.npz
echo ""

plotCorrelation -in multiBamSummary.results.npz --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap \
        --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCounts.png \
        --outFileCorMatrix SpearmanCorr_readCounts.tab
echo ""

######################################################
## STEP 5. plotCoverage: coverage of the sequencing ##
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
##################################################################################
## Step 6: FragmentSize :the fragment size distribution of your paired-end data ##
##################################################################################
echo "Step 6: fragment size and fingerprint"
if [ "${libraryType}" == "PE" ]; then
	for sample in ${sampleFiles[*]}
	do
		bamPEFragmentSize \
			-hist fragmentSize.${sample}.png \
			--plotTitle "Fragment Size" \
			--numberOfProcessors ${numberOfProcessors} \
			--maxFragmentLength 800 \
			--bamfiles ${sample}.dupMark.bam \
			--samplesLabel ${sample}
	echo ""
	done
fi

#plotFingerprint
plotFingerprint \
	--bamfiles ${dupMarkBamfiles} \
	--labels ${sampleFiles} \
	--minMappingQuality ${mappingQuality} \
	--skipZeros \
	--numberOfSamples 50000 \
	--numberOfProcessors ${numberOfProcessors} \
	--plotTitle "Fingerprints of samples" \
	--plotFile fingerprints.png \
	--outRawCounts fingerprints.tab \
	--outQualityMetrics qualityMetrics.tab
############################################################################
## STEP 7. coverage: to generate a sequencing-depth-normalized continuous ##
## profile of read coverages (BAM --> bigWig)                             ##
############################################################################
echo "Step 7: coverage."
for sample in ${sampleFiles[*]}
	do
	if [ "${libraryType}" == "SE" ]; then
                echo "generating normalized (RPKM) bigwig file of ${sample} w/ duplicate"
                bamCoverage --bam ${sample}.dupMark.bam --outFileName ${sample}.dupMark.rpkm.bw \
                        --binSize ${bwBinSize} \
                        --numberOfProcessors ${numberOfProcessors} --normalizeUsingRPKM \
                        --ignoreForNormalization chrM chrX chrY
                echo ""

                echo "generating unnormalized bigwig file of ${sample} w/ duplicate"
                bamCoverage --bam ${sample}.dupMark.bam --outFileName ${sample}.dupMark.bw \
                        --binSize ${bwBinSize} \
                        --numberOfProcessors ${numberOfProcessors}
                echo ""

                echo "generating normalized (RPKM) bigwig file of ${sample} w/o duplicate"
                bamCoverage  --bam ${sample}.dupMark.bam --outFileName ${sample}.dedup.rpkm.bw \
                        --binSize ${bwBinSize} --ignoreDuplicates \
                        --numberOfProcessors ${numberOfProcessors} --normalizeUsingRPKM \
                        --ignoreForNormalization chrM chrX chrY
                echo ""


                echo "generating unnormalized bigwig file of ${sample} w/o duplicate"
                bamCoverage --bam ${sample}.dupMark.bam --outFileName ${sample}.dedup.bw \
                        --binSize ${bwBinSize} --ignoreDuplicates \
                        --numberOfProcessors ${numberOfProcessors}
                echo ""
	elif [ "${libraryType}" == "PE" ]; then
		echo "generating normalized (RPKM) bigwig file of ${sample} w/ duplicate"
        	bamCoverage --bam ${sample}.dupMark.bam --outFileName ${sample}.dupMark.rpkm.bw \
                	--binSize ${bwBinSize} --extendReads \
                	--numberOfProcessors ${numberOfProcessors} --normalizeUsingRPKM \
                	--ignoreForNormalization chrM chrX chrY
        	echo ""

        	echo "generating unnormalized bigwig file of ${sample} w/ duplicate"
        	bamCoverage --bam ${sample}.dupMark.bam --outFileName ${sample}.dupMark.bw \
                	--binSize ${bwBinSize} --extendReads \
                	--numberOfProcessors ${numberOfProcessors}
        	echo ""

        	echo "generating normalized (RPKM) bigwig file of ${sample} w/o duplicate"
        	bamCoverage  --bam ${sample}.dedup.bam --outFileName ${sample}.dedup.rpkm.bw \
                	--binSize ${bwBinSize} --ignoreDuplicates --extendReads \
                	--numberOfProcessors ${numberOfProcessors} --normalizeUsingRPKM \
                	--ignoreForNormalization chrM chrX chrY
       		echo ""


        	echo "generating unnormalized bigwig file of ${sample} w/o duplicate"
        	bamCoverage --bam ${sample}.dedup.bam --outFileName ${sample}.dedup.bw \
                	--binSize ${bwBinSize} --ignoreDuplicates --extendReads \
                	--numberOfProcessors ${numberOfProcessors}
        	echo ""
	else
                echo "wrong library type"
	fi #fi is the concluding factor for bash if statements
done
