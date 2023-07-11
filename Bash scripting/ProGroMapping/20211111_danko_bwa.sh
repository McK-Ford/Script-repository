#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=groseq_%j.out
#SBATCH --mail-type=all
module load bwa
module load cutadapt
module load samtools
module load bedtools
module load seqtk
module load fastx-toolkit
module load prinseq
module load bedops
module load perl/5.10.1/b1

#########################################################################################
## About
## Preprocesses and aligns PRO-seq data. Takes PREFIX.fastq.gz (SE),  PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz (PE) or *.fastq.gz in the current working directory as input and writes BAM and bigWig files as output to the user-assigned output-dir.
## The output bigWig files ending with _minus.bw or _plus.bw are raw read counts without normalization. The RPM normalized outputs end with a suffix of .rpm.bw.
## Required options:
## --SEQ= SE or PE
## -i, --bwa-index=PATH   Path to the BWA index of the target genome
## -c, --chrom-info=PATH  Location of the chromInfo table.
## I/O options:
## -I, --fastq=PREFIX     Prefix for input files. Paired-end files require identical prefix and end with _R1.fastq.gz and _R2.fastq.gz"
## -T, --tmp=PATH         Path to a temporary storage directory.
## -O, --output-dir=DIR   Specify a directory to store output in.
## Required options for SE
## -G, --SE_READ=RNA_5prime Single-end sequencing from 5' end of nascent RNA, like GRO-seq.
## -P, --SE_READ=RNA_3prime Single-end sequencing from 3' end of nascent RNA, like PRO-seq.
## Options for PE
## --RNA5=R1_5prime    Specify the location of the 5' end of RNA
## --RNA3=R2_5prime    Specify the location of the 3' end of RNA"
## -5, --map5=TRUE     Report the 5' end of RNA [default on, --map5=TRUE].
## -3, --map5=FALSE    Report the 3' end of RNA, only available for PE [default off, --map5=TRUE].
## -s, --opposite-strand=TRUE
## --ADAPT_SE=TGGAATTCTCGGGTGCCAAGG 3' adapter to be removed from the 3' end of SE reads.       [default:TGGAATTCTCGGGTGCCAAGG]
## --ADAPT1=GATCGTCGGACTGTAGAACTCTGAACG 3' adapter to be removed from the 3' end of R2.      [default:GATCGTCGGACTGTAGAACTCTGAACG]
## --ADAPT2=AGATCGGAAGAGCACACGTCTGAACTC 3' adapter to be removed from the 3' end of R1. [default:AGATCGGAAGAGCACACGTCTGAACTC]
## --UMI1=0 The length of UMI barcode on the 5' of R1 read. [default: 0]
## --UMI2=0 The length of UMI barcode on the 5' of R2 read. [default: 0]
## When UMI1 or UMI2 are set > 0, the pipeline will perform PCR deduplicate.
## --Force_deduplicate=FALSE When --Force_deduplicate=TRUE, it will force the pipeline to       perform PCR deduplicate even there is no UMI barcode (i.e. UMI1=0 and UMI2=0). [default: FALSE]
## --ADD_B1=0         The length of additional barcode that will be trimmed on the 5' of R1 read. [default: 0]
## --ADD_B2=0         The length of additional barcode that will be trimmed on the 5' of R2 read. [default: 0]
## aln or mem algorithm?

files=("A_Rpn5" "A_Rpn9" "B_Rpn6" "C_Rpn7" "D_Rpn8")
SEQ="PE"
ADAPT_SE="TGGAATTCTCGGGTGCCAAGG"
ADAPT2="AGATCGGAAGAGCACACGTCTGAACTC"
ADAPT1="GATCGTCGGACTGTAGAACTCTGAACG"
UMI1=0
UMI2=0
Force_deduplicate=FALSE
ADD_B1=0
ADD_B2=0
map_L=0
OPP=FALSE
MAP5=TRUE
aln_mem="mem"
BWAIDX="/scratch/mford5/references/hg38/hg38.fa"
CHINFO="/scratch/mford5/references/hg38/hg38.chrom.sizes"
tmp="tmp"
mkdir $tmp
output="proseq_out"
mkdir $output
## INPUT & Parameters

if [[ "${files}" == "*.fastq.gz" ]]; then
    files=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| cut -d _ -f 2- |rev | sort | uniq`
fi

if [[ "$SEQ" == "PE" ]]; then
  if [[ -z "$RNA5" && -z "$RNA3" ]]; then
    RNA5="R1_5prime"
  fi
  if [[ "$RNA3" == "R1_5prime" ]]; then
   RNA5="R2_5prime"
  elif [[ "$RNA3" == "R2_5prime" ]]; then
    RNA5="R1_5prime"
  fi
fi

## Print out
echo "Processing PRO-seq data ..."
echo "SEQ mode is " $SEQ
echo "Files are " ${files[*]}
echo "Bwa index is " $BWAIDX
echo "chromInfo is " $CHINFO

if [[ "$SEQ" == "PE" ]] ; then
	## preprocesses data, removes adapters, trim.
	echo "Preprocessing fastq files:"
	mkdir ${tmp}/noadapt
	mkdir ${tmp}/passQC
	#do not currently use UMIs, need to add UMIs to this script.
	#there's also a function for getting num reads I need to look over.
	cutadapt -a ${ADAPT2} -e 0.10 --overlap 2 --output=${tmp}/${name}_trim_R1.fastq ${name}_R1.fastq.gz
	cutadapt -a ${ADAPT1} -e 0.10 --overlap 2 --output=${tmp}/${name}_trim_R2.fastq ${name}_R2.fastq.gz
	#Read 1, remove UMI2 and ADD_B2 from the 3' end of R1
	n2=$[UMI2+ADD_B2]
	cutadapt --cut -${n2} --minimum-length=10 ${tmp}/${name}_trim_R1.fastq --output=${tmp}/${name}_trim.${n2}Nremoved_R1.fastq -q 20 #remove the first n2 bases of the read, remove all reads shorter than 10, remove low quality ends
	cutadapt --minimum-length=10 ${name}_R1.fastq.gz --output=${tmp}/${name}_q20trim_R1.fastq -q 20
	cat ${tmp}/${name}_q20trim_R1.fastq ${tmp}/${name}_trim.${n2}Nremoved_R1.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${tmp} -k1,1 -S 10G | tr '\t' '\n' > ${tmp}/noadapt/${name}_noadapt_R1.fastq
	# Read2
	# remove UMI1 and ADD_B1 from the 3 prime end of R2
	n1=$[UMI1+ADD_B1]
	cutadapt --cut -${n1} --minimum-length=10 ${tmp}/${name}_trim_R2.fastq --output=${tmp}/${name}_trim.${n1}Nremoved_R2.fastq -q 20
	cutadapt --minimum-length=10 ${name}_R2.fastq.gz --output=${tmp}/${name}_q20trim_R2.fastq -q 20
	cat ${tmp}/${name}_q20trim_R2.fastq ${tmp}/${name}_trim.${n1}Nremoved_R2.fastq | paste - - - - | LC_ALL=C sort --temporary-directory=${tmp} -k1,1 -S 10G | tr '\t' '\n' > ${tmp}/noadapt/${name}_noadapt_R2.fastq
	#get counts
	echo 'Number of reads after adapter removal and QC:' >> ${tmp}/${name}.QC.log
	echo "R1" >> ${tmp}/${name}.QC.log
	cat ${tmp}/noadapt/${name}_noadapt_R1.fastq | grep @ -c >> ${tmp}/${name}.QC.log
	echo "R2" >> ${tmp}/${name}.QC.log
	cat ${tmp}/noadapt/${name}_noadapt_R2.fastq | grep @ -c >> ${tmp}/${name}.QC.log
	rm ${tmp}/${name}_trim_R1.fastq ${tmp}/${name}_q20trim_R1.fastq ${tmp}/${name}_trim.${n2}Nremoved_R1.fastq
	rm ${tmp}/${name}_trim_R2.fastq ${tmp}/${name}_q20trim_R2.fastq ${tmp}/${name}_trim.${n1}Nremoved_R2.fastq
	## Collapse reads using prinseq-lite.pl. if there are UMI barcodes or $Force_deduplicate=TRUE
	    if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then # if there is UMI barcode, Will perform deduplicates with first ${dedup_L} bp
	      # Remove PCR duplciates.
	      # fastq with the first dedup_L nt
	      cat ${tmp}/noadapt/${name}_noadapt_R1.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq
	      cat ${tmp}/noadapt/${name}_noadapt_R2.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq
	       # deduplicate using the first dedup_L nt
	       prinseq-lite.pl -fastq ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq -fastq2 ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq -derep 1 -out_format 3 -out_bad null -out_good ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup -min_len 15 2> ${output}/${name}.prinseq-pcrDups.gd
	       rm ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R1.fastq ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}_R2.fastq
	       # make a list of name from deduplicated fastq
	       cat ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1.fastq | awk '(NR%4==1){print substr($1,2)}' > ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt
	       cat ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons.fastq | awk '(NR%4==1){print substr($1,2)}' > ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons_l${dedup_L}.txt
	       cat ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons.fastq | awk '(NR%4==1){print substr($1,2)}' > ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons_l${dedup_L}.txt
	       rm ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons.fastq
	       # generate fastq from the list of name
	       seqtk subseq ${tmp}/noadapt/${name}_noadapt_R1.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${tmp}/passQC/${name}_dedup_withBarcode_1.fastq
	       seqtk subseq ${tmp}/noadapt/${name}_noadapt_R2.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${tmp}/passQC/${name}_dedup_withBarcode_2.fastq
	       seqtk subseq ${tmp}/noadapt/${name}_noadapt_R1.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_1_singletons_l${dedup_L}.txt > ${tmp}/passQC/${name}_dedup_1_singletons.fastq
	       seqtk subseq ${tmp}/noadapt/${name}_noadapt_R2.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_2_singletons_l${dedup_L}.txt > ${tmp}/passQC/${name}_dedup_2_singletons.fastq
	       rm ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup*l${dedup_L}.txt
	       # trim the UMI and additional barcode after dereplicate
	       prinseq-lite.pl -trim_left ${n1} -fastq ${tmp}/passQC/${name}_dedup_withBarcode_1.fastq -  out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1 2>> ${output}/${name}.prinseq-pcrDups.gd
	       prinseq-lite.pl -trim_left ${n2} -fastq ${tmp}/passQC/${name}_dedup_withBarcode_2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2 2>> ${output}/${name}.prinseq-pcrDups.gd
	       rm ${tmp}/passQC/${name}_dedup_withBarcode_1.fastq ${tmp}/passQC/${name}_dedup_withBarcode_2.fastq
	       # min_len 15
	       prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1.fastq -fastq2 ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_QC_end 2>> ${output}/${name}.prinseq-pcrDups.gd
	       rm ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1.fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2.fastq
	       echo 'Number of paired reads after PCR duplicates removal and QC:' >> ${tmp}/${name}.QC.log
	       cat ${tmp}/passQC/${name}_dedup_QC_end_1.fastq | grep @ -c >> ${tmp}/${name}.QC.log
	   elif [[ ${Force_deduplicate} == "TRUE" ]]; then # if there is NO UMI barcode, Will perform deduplicates with whole legnth of reads
	      # Remove PCR duplciates.
	      prinseq-lite.pl -derep 1 -fastq ${tmp}/noadapt/${name}_noadapt_R1.fastq -fastq2 ${tmp}/noadapt/${name}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_withBarcode 2> ${output}/${name}.prinseq-pcrDups.gd
	      # trim the UMI and additional barcode after dereplicate
	      prinseq-lite.pl -trim_left ${n1} -fastq ${tmp}/passQC/${name}_dedup_withBarcode_1.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1 2>> ${output}/${name}.prinseq-pcrDups.gd
	      prinseq-lite.pl -trim_left ${n2} -fastq ${tmp}/passQC/${name}_dedup_withBarcode_2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2 2>> ${output}/${name}.prinseq-pcrDups.gd
	      rm ${tmp}/passQC/${name}_dedup_withBarcode_1.fastq ${tmp}/passQC/${name}_dedup_withBarcode_2.fastq
	      # min_len 15
	      prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1.fastq -fastq2 ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_QC_end 2>> ${output}/${name}.prinseq-pcrDups.gd
	      rm ${tmp}/passQC/${name}_dedup_BarcodeRemoved_1.fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved_2.fastq
	      echo 'Number of paired reads after PCR duplicates removal and QC:' >> ${tmp}/${name}.QC.log
	      cat ${tmp}/passQC/${name}_dedup_QC_end_1.fastq | grep @ -c >> ${tmp}/${name}.QC.log
	  else #THIS HAS BEEN TESTED CAUSE I DIDN'T LIST UMIs
	     # trim the additional barcode ${ADD_B1} and ${ADD_B2} WITHOUT dereplicate. If no barcode, prinseq-lite.pl will remove unpair reads and reads that are length 0
	     prinseq-lite.pl -trim_left ${ADD_B1} -fastq ${tmp}/noadapt/${name}_noadapt_R1.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_BarcodeRemoved_1 2>> ${output}/${name}.prinseq-pcrDups.gd #trim 5' by add_b1 if desired, get fastq as output.
	     prinseq-lite.pl -trim_left ${ADD_B2} -fastq ${tmp}/noadapt/${name}_noadapt_R2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_BarcodeRemoved_2 2>> ${output}/${name}.prinseq-pcrDups.gd
	     # min_len 15
	     prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_BarcodeRemoved_1.fastq -fastq2 ${tmp}/passQC/${name}_BarcodeRemoved_2.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_QC_end 2> ${output}/${name}.prinseq-pcrDups.gd
	     rm ${tmp}/passQC/${name}_BarcodeRemoved_1.fastq ${tmp}/passQC/${name}_BarcodeRemoved_2.fastq
	     echo 'Number of paired reads after final QC:' >> ${tmp}/${name}.QC.log
	     cat ${tmp}/passQC/${name}_QC_end_1.fastq | grep @ -c >> ${tmp}/${name}.QC.log
	 fi
	 cat ${output}/${name}.prinseq-pcrDups.gd
	done
	# if there is a data-set wide length cutoff map_L
	if [[ ${map_L} != 0 ]]; then
		 mkdir ${tmp}/passQC_length_${map_L}
		 for f in ${tmp}/passQC/*_QC_end_1.fastq
			do name=`echo $f |rev |cut -d / -f 1 |cut -d _ -f 4-|rev`
			cat ${f} | fastx_trimmer -Q33 -l ${map_L} -o ${tmp}/passQC_length_${map_L}/${name}_l${map_L}_QC_end_1.fastq
		 done
		 for f in ${tmp}/passQC/*_QC_end_2.fastq
			do name=`echo $f |rev |cut -d / -f 1 |cut -d _ -f 4-|rev`
			cat ${f} | fastx_trimmer -Q33 -l ${map_L} -o ${tmp}/passQC_length_${map_L}/${name}_l${map_L}_QC_end_2.fastq
		 done
	fi
	fi
	## Cleanup.
	# rm -r ${tmp}/noadapt
	gzip ${tmp}/passQC*/*.fastq
#############################################
if [[ ${map_L} != 0 ]]; then
   ToMapDir=${tmp}/passQC_length_${map_L}
else
   ToMapDir=${tmp}/passQC
fi
echo "Reads for mapping is from : $ToMapDir"
#QC_INPUT=ls  ${ToMapDir}/*_QC_end_1.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq
# main map
echo "Mapping reads:"

if [${aln_mem}=="aln"]; then # elif -aln was given, use aln
 for name in ${files[*]}
  do
  ## Align using BWA aln
  bwa aln ${BWAIDX} ${ToMapDir}/${name}_R1.fastq.gz  >  ${tmp}/${name}_aln_sa1.sai
  bwa aln ${BWAIDX} ${ToMapDir}/${name}_R2.fastq.gz  >  ${tmp}/${name}_aln_sa2.sai
  bwa sampe -n 1 -f ${ToMapDir}/${name}_end.sam ${BWAIDX} ${tmp}/${name}_aln_sa1.sai ${tmp}/${name}_aln_sa2.sai ${ToMapDir}/${name}_R1.fastq.gz ${ToMapDir}/${name}_R2.fastq.gz
  samtools view -bf 0x2 -q 20 ${ToMapDir}/${name}_end.sam | samtools sort -n - > ${tmp}/${name}.sort.bam
  rm ${tmp}/${name}_aln_sa1.sai ${tmp}/${name}_aln_sa2.sai ${ToMapDir}/${name}_end.sam
  done
else #defaul use mem
  for name in ${files[*]}
    do
    ## Align using BWA.
	echo "using mem"
    bwa mem -k 19 ${BWAIDX} ${ToMapDir}/${name}_R1.fastq.gz ${ToMapDir}/${name}_R2.fastq.gz | \
    samtools view -bf 0x2 -q 20 - | samtools sort -n - > ${tmp}/${name}.sort.bam
  done
fi
for name in ${files[*]}
  do
  cp ${tmp}/${name}.sort.bam ${output}
done
## Cleanup
#find ${tmp} -name "*.sort.bam" -size -1024k -delete
##########################################
## Write out the bigWigs.
echo "Writing bigWigs:"
for bam in ${tmp}/*.sort.bam
 do
   newname=`echo $bam | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev`
   echo ${newname} > ${output}/${newname}.align.log
	bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' | gzip > ${tmp}/${newname}.RNA5.bed.gz
	bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | gzip > ${tmp}/${newname}.long.bedpe.gz
zcat ${tmp}/${newname}.regionsfilt.bedpe.gz | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$2,$6,$7,$8,$9}; ($10 == "-") {print $1,$2,$6,$7,$8,$9}'|gzip > ${tmp}/${newname}.long.bed.gz
zcat ${tmp}/${newname}.long.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.long.bed.gz
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.long.bed.gz -g ${CHINFO} -strand + > ${tmp}/${newname}_long_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.long.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}_long_minus.bedgraph
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/${newname}_long_plus.bedgraph ${CHINFO} ${tmp}/${newname}_long_plus.bw
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/${newname}_long_minus.bedgraph ${CHINFO} ${tmp}/${newname}_long_minus.bw

	bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}'|gzip > ${tmp}/${newname}.RNA3.bed.gz
readCount='zcat ${tmp}/${genename}.RNA5.bed.gz'
#remove rRNA and reverse the strand
zcat ${tmp}/${newname}.RNA5.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.RNA5.bed.gz
zcat ${tmp}/${newname}.RNA3.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.RNA3.bed.gz
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.RNA5.bed.gz -g ${CHINFO} -strand + > ${tmp}/${newname}_RNA5_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.RNA5.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}_RNA5_minus.bedgraph
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.RNA3.bed.gz -g ${CHINFO} -strand + > ${tmp}/${newname}_RNA3_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.RNA3.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}_RNA3_minus.bedgraph
cat ${tmp}/${newname}_RNA5_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_RNA5_plus.rpm.bedGraph}
cat ${tmp}/${newname}_RNA5_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_RNA5_minus.rpm.bedGraph}
cat ${tmp}/${newname}_RNA3_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_RNA3_plus.rpm.bedGraph}
cat ${tmp}/${newname}_RNA3_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_RNA3_minus.rpm.bedGraph}
done
for bam in ${tmp}/*.sort.bam
 do
 newname=`echo $bam | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev`
echo ${newname} > ${output}/${newname}.align.log
bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' | gzip > ${tmp}/${newname}.RNA5.bed.gz


	bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | gzip > ${tmp}/${newname}.long.bedpe.gz

zcat ${tmp}/${newname}.regionsfilt.bedpe.gz | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}'|gzip > ${tmp}/${newname}.TSS.bed.gz
zcat ${tmp}/${newname}.regionsfilt.bedpe.gz | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$6,$6+1,$7,$8,$10}; ($10 == "-") {print $1,$2,$2+1,$7,$8,$10}'|gzip > ${tmp}/${newname}.Pauseloc.bed.gz
zcat ${tmp}/${newname}.long.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.long.bed.gz
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.long.bed.gz -g ${CHINFO} -strand + > ${tmp}/${newname}_long_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.long.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}_long_minus.bedgraph
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/${newname}_long_plus.bedgraph ${CHINFO} ${tmp}/${newname}_long_plus.bw
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/${newname}_long_minus.bedgraph ${CHINFO} ${tmp}/${newname}_long_minus.bw

bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}'|gzip > ${tmp}/${newname}.RNA3.bed.gz
done
zcat ${tmp}/*.TSS.bed.gz | cat | gzip > ${tmp}/TSS.bed.gz
zcat ${tmp}/*.Pauseloc.bed.gz | cat | gzip > ${tmp}/Pauseloc.bed.gz
zcat ${tmp}/TSS.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/TSS.nr.rs.bed.gz
zcat ${tmp}/Pauseloc.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/Pauseloc.nr.rs.bed.gz
bedtools genomecov -bg -i ${tmp}/TSS.nr.rs.bed.gz -g ${CHINFO} -strand + > ${tmp}/TSS_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/TSS.nr.rs.bed.gz -g ${CHINFO} -strand - > ${tmp}/TSS_minus.bedgraph
bedtools genomecov -bg -i ${tmp}/Pauseloc.nr.rs.bed.gz -g ${CHINFO} -strand + > ${tmp}/Pauseloc_plus.bedgraph
bedtools genomecov -bg -i ${tmp}/Pauseloc.nr.rs.bed.gz -g ${CHINFO} -strand - > ${tmp}/Pauseloc_minus.bedgraph
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/Pauseloc_plus.bedgraph ${CHINFO} ${tmp}/Pauseloc_plus.bw
/scratch/mford5/tools/bedGraphToBigWig ${tmp}/Pauseloc_minus.bedgraph ${CHINFO} ${tmp}/Pauseloc_minus.bw
fi #*/
