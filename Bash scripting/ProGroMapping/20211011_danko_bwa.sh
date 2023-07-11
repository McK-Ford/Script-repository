#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=1
#SBATCH -output=groseq_%j.out
#SBATCH --mail-type=all
module load bwa
module load cutadapt
module load samtools
module load seqtk
module load fastx_toolkit
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

files = (A_Rpn5 A_Rpn9 B_Rpn6 C_Rpn7 D_Rpn8)
SEQ = "PE"
ADAPT_SE ="TGGAATTCTCGGGTGCCAAGG"
ADAPT2 ="AGATCGGAAGAGCACACGTCTGAACTC"
ADAPT1 ="GATCGTCGGACTGTAGAACTCTGAACG"
UMI1 = 0
UMI2 = 0
Force_deduplicate == FALSE
ADD_B1 = 0
ADD_B2 = 0
map_L = 0
OPP = FALSE
MAP5 = TRUE
aln_mem = "mem"
bwaIndex = /scratch/mford5/references/hg38/hg38.fa
chromInfo = /scratch/mford5/references/hg38/hg38.chrom.sizes
tmp = "tmp"
mkdir $tmp
output = "proseq_out"
mkdir $output
## INPUT & Parameters
# PE
if [[ "$files" == "*.fastq.gz" ]]; then
    files=`ls *.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3-| cut -d _ -f 2- |rev | sort | uniq`
fi
# SE
if [[ "$SEQ" == "SE" ]] ; then
  if [[ "$SE_READ" == "RNA_5prime" ]] ; then
    SE_output="G"
    RNA5="R1_5prime"
    OPP="FALSE"
  elif [[ "$SE_READ" == "RNA_3prime" ]] ; then
    SE_output="P"
    RNA3="R1_5prime"
    OPP="TRUE"
  fi
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
echo "Files are " $files
echo "Bwa index is " $BWAIDX
echo "chromInfo is " $CHINFO

if [[ "$SEQ" == "PE" ]] ; then
  #############################################
  ## Preprocess data.  Remove adapters.  Trim.
  echo "Preprocessing fastq files:"
  mkdir ${tmp}/noadapt
  mkdir ${tmp}/passQC
  if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then
    dedup_L=30
    mkdir ${tmp}/noadapt/l${dedup_L}
    mkdir ${tmp}/noadapt/l${dedup_L}_nodups
  fi
  for name in ${files[*]}
  do
    ## get num reads
    echo 'Number of original input reads:' > ${tmp}/${name}.QC.log
    echo ${name}_R1.fastq.gz >> ${tmp}/${name}.QC.log
    zcat ${name}_R1.fastq.gz | grep @ -c >> ${tmp}/${name}.QC.log
    echo ${name}_R2.fastq.gz >> ${tmp}/${name}.QC.log
    zcat ${name}_R2.fastq.gz | grep @ -c >> ${tmp}/${name}.QC.log
    ## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
    # Remove adapter
    cutadapt -a ${ADAPT2} -e 0.10 --overlap 2 --output=${tmp}/${name}_trim_R1.fastq ${name}_R1.fastq.gz #trim adaptor with max error rate of .10, need at least two bases of the adaptor to trim,
    cutadapt -a ${ADAPT1} -e 0.10 --overlap 2 --output=${tmp}/${name}_trim_R2.fastq ${name}_R2.fastq.gz
    # Read1
    # remove UMI2 and ADD_B2 from the 3 prime end of R1
    n2=$[UMI2+ADD_B2]
    cutadapt --cut -${n2} --minimum-length=10 ${tmp}/${name}_trim_R1.fastq --output=${tmp}/${name}_trim.${n2}Nremoved_R1.fastq -q 20 #remove the first n2 bases of the read, remove all reads shorter than 10, remove low quality ends
    cutadapt --minimum-length=10 ${tmp}/${name}_untrim_R1.fastq --output=${tmp}/${name}_q20trim_R1.fastq -q 20
    cat ${tmp}/${name}_q20trim_R1.fastq ${tmp}/${name}_trim.${n2}Nremoved_R1.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${tmp} -k1,1 -S 10G | tr '\t' '\n' > ${tmp}/noadapt/${name}_noadapt_R1.fastq
    # Read2
    # remove UMI1 and ADD_B1 from the 3 prime end of R2
    n1=$[UMI1+ADD_B1]
    cutadapt --cut -${n1} --minimum-length=10 ${tmp}/${name}_trim_R2.fastq --output=${tmp}/${name}_trim.${n1}Nremoved_R2.fastq -q 20
    cutadapt --minimum-length=10 ${tmp}/${name}_untrim_R2.fastq --output=${tmp}/ ${name}_q20trim_R2.fastq -q 20
    cat ${tmp}/${name}_q20trim_R2.fastq ${tmp}/${name}_trim.${n1}Nremoved_R2.fastq | paste - - - - | LC_ALL=C sort --temporary-directory=${tmp} -k1,1 -S 10G | tr '\t' '\n' > ${tmp}/noadapt/${name}_noadapt_R2.fastq
    #get counts
    echo 'Number of reads after adapter removal and QC:' >> ${tmp}/${name}.QC.log
    echo "R1" >> ${tmp}/${name}.QC.log
    cat ${tmp}/noadapt/${name}_noadapt_R1.fastq | grep @ -c >> ${tmp}/${name}.QC.log
    echo "R2" >> ${tmp}/${name}.QC.log
    cat ${tmp}/noadapt/${name}_noadapt_R2.fastq | grep @ -c >> ${tmp}/${name}.QC.log
    rm ${tmp}/${name}_trim_R1.fastq ${tmp}/${name}_untrim_R1.fastq ${tmp}/${name}_q20trim_R1.fastq ${tmp}/${name}_trim.${n2}Nremoved_R1.fastq
    rm ${tmp}/${name}_trim_R2.fastq ${tmp}/${name}_untrim_R2.fastq ${tmp}/${name}_q20trim_R2.fastq ${tmp}/${name}_trim.${n1}Nremoved_R2.fastq
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
  else
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
## Cleanup.
rm -r ${tmp}/noadapt
gzip ${tmp}/passQC*/*.fastq

#############################################
## Align reads.
if [[ ${map_L} != 0 ]]; then
   ToMapDir=${tmp}/passQC_length_${map_L}
else
   ToMapDir = ${tmp}/passQC
fi
echo "Reads for mapping is from : $ToMapDir"
QC_INPUT=`ls  ${ToMapDir}/*_QC_end_1.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq`
# main map
echo "Mapping reads:"

if [aln_mem == "aln"]; then # elif -aln was given, use aln
 for name in ${QC_INPUT}
  do
  ## Align using BWA aln
  bwa aln ${BWAIDX} ${ToMapDir}/${name}_1.fastq.gz  >  ${tmp}/${name}_aln_sa1.sai
  bwa aln ${BWAIDX} ${ToMapDir}/${name}_2.fastq.gz  >  ${tmp}/${name}_aln_sa2.sai
  bwa sampe -n 1 -f ${ToMapDir}/${name}_end.sam ${BWAIDX} ${tmp}/${name}_aln_sa1.sai ${tmp}/${name}_aln_sa2.sai ${ToMapDir}/${name}_1.fastq.gz ${ToMapDir}/${name}_2.fastq.gz
  samtools view -bf 0x2 -q 20 ${ToMapDir}/${name}_end.sam | samtools sort -n - > ${tmp}/${name}.sort.bam
  rm ${tmp}/${name}_aln_sa1.sai ${tmp}/${name}_aln_sa2.sai ${ToMapDir}/${name}_end.sam
  done
else #defaul use mem
  for name in ${QC_INPUT}
    do
    ## Align using BWA.
    bwa mem -k 19 ${BWAIDX} ${ToMapDir}/${name}_1.fastq.gz ${ToMapDir}/${name}_2.fastq.gz | \
    samtools view -bf 0x2 -q 20 - | samtools sort -n - > ${tmp}/${name}.sort.bam
  done
fi
for name in ${QC_INPUT}
  do
  cp ${tmp}/${name}.sort.bam ${output}
done
## Cleanup
find ${tmp} -name "*.sort.bam" -size -1024k -delete
#############################################
## Write out the bigWigs.
echo "Writing bigWigs:"
for bam in ${tmp}/*.sort.bam
 do
   newname =`echo $bam | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
   echo ${newname} > ${output}/${newname}.align.log
   if [ "${RNA5}" == "R1_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA. Danko lab leChRO-Seq protocol is on the 5' of _R1 readl, same strand of R1 ($9)
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' | gzip > ${tmp}/${newname}.bed.gz
       else ## report The 3' end of the RNA.  Danko lab leChRO-Seq protocol is on the 5 prime of _R2 read, opposite strand of R2 (R2 strand $10, R1 strand $9)
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}' | gzip > ${tmp}/${newname}.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then  ## report The 5' end of the RNA.
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' | gzip > ${tmp}/${newname}.bed.gz
       else ## report The 3' end of the RNA.
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}; ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}' | gzip > ${tmp}/${newname}.bed.gz
       fi
     fi
   elif [ "${RNA5}" == "R2_5prime" ] ; then
     if [ "${OPP}" == "FALSE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$10}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$10}'|gzip > ${tmp}/${newname}.bed.gz
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$10}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$10}' |gzip  > ${tmp}/${newname}.bed.gz
       fi
     elif [ "${OPP}" == "TRUE" ] ; then
       if [ "${MAP5}" == "TRUE" ] ; then #report the 5 prime end of RNA, in Engreitz data is 5 prime end of R2, same strand
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/$${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($10 == "+") {print $1,$5,$5+1,$7,$8,$9}; ($10 == "-") {print $1,$6-1,$6,$7,$8,$9}'|gzip > ${tmp}/${newname}.bed.gz
       else  ## report the 3-prime end of the RNA, in Engreitz data is the 5' end of R1 read, but opposite strand
          bedtools bamtobed -bedpe -mate1 -i $bam 2> ${tmp}/${newname}.kill.warnings | awk 'BEGIN{OFS="\t"} ($9 == "+") {print $1,$2,$2+1,$7,$8,$9}; ($9 == "-") {print $1,$3-1,$3,$7,$8,$9}' |gzip  > ${tmp}/${newname}.bed.gz
       fi
     fi
   fi

   echo 'Number of mappable reads:' >> ${output}/${newname}.align.log
   readCount=`zcat ${tmp}/${newname}.bed.gz | grep "" -c`
   echo ${readCount} >> ${output}/${newname}.align.log
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${tmp}/${newname}.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${output}/${newname}.align.log
   echo `zcat ${tmp}/${newname}.nr.rs.bed.gz | grep "" -c` >> ${output}/${newname}.align.log
   ## Convert to bedGraph ... Can't gzip these, unfortunately.
   bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.bed.gz -g ${CHINFO} -strand + > ${tmp}/$j\_plus.bedGraph
   bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}_minus.noinv.bedGraph
   ## Invert minus strand.
   cat ${tmp}/${newname}_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${tmp}/${newname}_minus.bedGraph ## Invert read counts on the minus strand.
   ## normalized by RPM
   cat ${tmp}/${newname}_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_plus.rpm.bedGraph
   cat ${tmp}/${newname}_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_minus.rpm.bedGraph
   ## Then to bigWig (nomalized and non-nomrmalized ones)
   bedGraphToBigWig ${tmp}/${newname}_plus.rpm.bedGraph ${CHINFO} ${output}/${newname}_plus.rpm.bw
   bedGraphToBigWig ${tmp}/${newname}_minus.rpm.bedGraph ${CHINFO} ${output}/${newname}_minus.rpm.bw
   bedGraphToBigWig ${tmp}/${newname}_plus.bedGraph ${CHINFO} ${output}/${newname}_plus.bw
   bedGraphToBigWig ${tmp}/${newname}_minus.bedGraph ${CHINFO} ${output}/${newname}_minus.bw
   rm ${tmp}/${newname}.nr.rs.bed.gz ${tmp}/${newname}.bed.gz ${tmp}/${newname}*.bedGraph
 done
fi #*/

if [[ "$SEQ" == "SE" ]] ; then
#############################################
## Preprocess data.  Remove adapters.  Trim.
echo "Preprocessing fastq files:"
mkdir ${tmp}/noadapt
mkdir ${tmp}/passQC
if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then
  dedup_L=30
  mkdir ${tmp}/noadapt/l${dedup_L}
  mkdir ${tmp}/noadapt/l${dedup_L}_nodups
fi
for name in ${files[*]}
 do
  ## read stats
  echo ${name} >  ${tmp}/${name}.QC.log
  echo 'Number of original input reads:' >> ${tmp}/${name}.QC.log
  zcat ${name}.fastq.gz | grep @ -c >> ${tmp}/${name}.QC.log
  ## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
  # Remove adapter
  cutadapt -a ${ADAPT_SE} -e 0.10 --overlap 2 --output=${tmp}/${name}_trim.fastq --untrimmed-output=${tmp}/${name}_untrim.fastq ${old_name}.fastq.gz
  # Read1
  # remove UMI2 and ADD_B2 from the 3 prime end of R1
   n2=$[UMI2+ADD_B2]
   n1=$[UMI1+ADD_B1]
   cutadapt --cut -${n2} --minimum-length=10 ${tmp}/${name}_trim.fastq --output=${tmp}/${name}_trim.${n2}Nremoved.fastq -q 20
   cutadapt --minimum-length=10 ${tmp}/${name}_untrim.fastq --output=${tmp}/${name}_q20trim.fastq -q 20
   cat ${tmp}/${name}_q20trim.fastq ${tmp}/${name}_trim.${n2}Nremoved.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${tmp} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${tmp}/noadapt/${name}_noadapt.fastq

  echo 'Number of reads after adapter removal and QC:' >> ${tmp}/${name}.QC.log
  cat ${tmp}/noadapt/${name}_noadapt.fastq | grep @ -c >> ${tmp}/${name}.QC.log

  rm ${tmp}/${name}_trim.fastq ${tmp}/${name}_untrim.fastq ${tmp}/${name}_q20trim.fastq ${tmp}/${name}_trim.${n2}Nremoved.fastq

  ## Collapse reads using prinseq-lite.pl. if there are UMI barcodes
  if [[ ${UMI2} != 0 || ${UMI1} != 0 ]]; then
   # Remove PCR duplciates.
   # fastq with the first dedup_L nt
   cat ${tmp}/noadapt/${name}_noadapt.fastq | fastx_trimmer -Q33 -l ${dedup_L} -o ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq
   # deduplicate using the first dedup_L nt
   prinseq-lite.pl -fastq ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq -derep 1 -out_format 3 -out_bad null -out_good ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup -min_len 15 2> ${output}/${name}.prinseq-pcrDups.gd
   rm ${tmp}/noadapt/l${dedup_L}/${name}_noadapt_l${dedup_L}.fastq
   # make a list of name from deduplicated fastq
   cat ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup.fastq | awk '(NR%4==1){print substr($1,2)}' > ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt
   rm ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup.fastq
   # generate fastq from the list of name
   seqtk subseq ${tmp}/noadapt/${name}_noadapt.fastq ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt > ${tmp}/passQC/${name}_dedup_withBarcode.fastq
   rm ${tmp}/noadapt/l${dedup_L}_nodups/${name}_dedup_l${dedup_L}.txt
   # trim the UMI and additional barcode after dereplicate
   prinseq-lite.pl -trim_left ${n1} -fastq ${tmp}/passQC/${name}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved 2>> ${output}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_QC_end 2>> ${output}/${name}.prinseq-pcrDups.gd
   rm ${tmp}/passQC/${name}_dedup_withBarcode.fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved.fastq
   echo 'Number of reads after PCR duplicates removal and QC:' >> ${tmp}/${name}.QC.log
   cat ${tmp}/passQC/${name}_dedup_QC_end.fastq | grep @ -c >> ${tmp}/${name}.QC.log
  elif [[ ${Force_deduplicate} == "TRUE" ]]; then
   # Remove PCR duplciates.
   prinseq-lite.pl -derep 1 -fastq ${tmp}/noadapt/${name}_noadapt.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_withBarcode 2> ${output}/${name}.prinseq-pcrDups.gd
   # trim the UMI and additional barcode after dereplicate
   prinseq-lite.pl -trim_left ${n1} -fastq ${tmp}/passQC/${name}_dedup_withBarcode.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_BarcodeRemoved 2>> ${output}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_dedup_QC_end 2>> ${output}/${name}.prinseq-pcrDups.gd
   rm ${tmp}/passQC/${name}_dedup_withBarcode.fastq ${tmp}/passQC/${name}_dedup_BarcodeRemoved.fastq
   echo 'Number of reads after PCR duplicates removal and QC:' >> ${tmp}/${name}.QC.log
   cat ${tmp}/passQC/${name}_dedup_QC_end.fastq | grep @ -c >> ${tmp}/${name}.QC.log
  else
   # trim the additional barcode WITHOUT dereplicate. If no barcode, prinseq-lite.pl will remove unpair reads and reads that are length 0
   prinseq-lite.pl -trim_left ${ADD_B1} -fastq ${tmp}/noadapt/${name}_noadapt.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_BarcodeRemoved 2>> ${output}/${name}.prinseq-pcrDups.gd
   # min_len 15
   prinseq-lite.pl -min_len 15 -fastq ${tmp}/passQC/${name}_BarcodeRemoved.fastq -out_format 3 -out_bad null -out_good ${tmp}/passQC/${name}_QC_end 2> ${output}/${name}.prinseq-pcrDups.gd
   rm ${tmp}/passQC/${name}_BarcodeRemoved.fastq
   echo 'Number of reads after final QC:' >> ${tmp}/${name}.QC.log
   cat ${tmp}/passQC/${name}_QC_end.fastq | grep @ -c >> ${tmp}/${name}.QC.log
  fi
  cat ${output}/${name}.prinseq-pcrDups.gd
done

# if there is a data-set wide length cutoff map_L
if [[ ${map_L} != 0 ]]; then
   mkdir ${tmp}/passQC_length_${map_L}
   for f in ${tmp}/passQC/*_QC_end.fastq
    do name=`echo $f |rev |cut -d / -f 1 |cut -d _ -f 3-|rev`
    cat ${f} | fastx_trimmer -Q33 -l ${map_L} -o ${tmp}/passQC_length_${map_L}/${name}_l${map_L}_QC_end.fastq
   done
fi

## Cleanup.
rm -r ${tmp}/noadapt
gzip ${tmp}/passQC*/*.fastq

#############################################
## Align reads.

if [[ ${map_L} != 0 ]]; then
   ToMapDir=${tmp}/passQC_length_${map_L}
else
   ToMapDir=${tmp}/passQC
fi

QC_INPUT=`ls  ${ToMapDir}/*_QC_end.fastq.gz | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- | cut -d _ -f 2- | rev| sort | uniq`

echo " "
echo "Mapping reads:"
echo "Reads for mapping is from :  $ToMapDir"
# main map
  if [aln_mem == "mem"]; then #use mem
    for name in ${QC_INPUT}
      do
      ## Align using BWA.
      bwa mem -k 19 ${BWAIDX} ${ToMapDir}/${name}_end.fastq.gz | \
      samtools view -b -q 20 - | samtools sort -n - > ${tmp}/${name}.sort.bam
      done
    for name in ${QC_INPUT}
      do
      cp ${tmp}/${name}.sort.bam ${output}  ## Saves the sorted BAM in the output file.
      done
  else #default use aln
    for name in ${QC_INPUT}
      do
      ## Align using BWA aln
echo "bwa aln ${BWAIDX} ${ToMapDir}/${name}_end.fastq.gz | \
      bwa samse -n 1 -f ${ToMapDir}/${name}_end.sam ${BWAIDX} - ${ToMapDir}/${name}_end.fastq.gz "
      bwa aln ${BWAIDX} ${ToMapDir}/${name}_end.fastq.gz | \
      bwa samse -n 1 -f ${ToMapDir}/${name}_end.sam ${BWAIDX} - ${ToMapDir}/${name}_end.fastq.gz
    done
    for name in ${QC_INPUT}
      do
      echo "samtools view -b -q 20 ${ToMapDir}/${name}_end.sam | samtools sort -n - > ${tmp}/${name}.sort.bam"
      samtools view -b -q 20 ${ToMapDir}/${name}_end.sam | samtools sort -n - > ${tmp}/${name}.sort.bam
    done
    for name in ${QC_INPUT}
      do
      rm ${ToMapDir}/${name}_end.sam
      cp ${tmp}/${name}.sort.bam ${output}  ## Saves the sorted BAM in the output file.
    done
   fi
## Cleanup
find ${tmp} -name "*.sort.bam" -size -1024k -delete

#############################################
## Write out the bigWigs.
echo "Writing bigWigs:"
for bam in ${tmp}/*.sort.bam
do
   newname=`echo $bam | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
   echo $newname > ${output}/${newname}.align.log
# in SE, MAP5 alwasys TRUE
   #if [[ "${RNA5}" == "R1_5prime" && "${OPP}" == "FALSE" ]] ; then ## report The 5 prime end of the RNA.   #like GRO-seq
   if [[ "$SE_output" == "G" ]] ; then
      bedtools bamtobed -i $bam 2> ${tmp}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > ${tmp}/${newname}.bed.gz
   #elif [[ "${RNA3}" == "R1_5prime" && "${OPP}" == "TRUE" ]] ; then  #like PRO-seq
    elif [[ "$SE_output" == "P" ]] ; then
      bedtools bamtobed -i $bam 2> ${tmp}/kill.warnings| awk 'BEGIN{OFS="\t"} ($5 > 0){print $0}' | \
      awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,"-"}; ($6 == "-") {print $1,$3-1,$3,$4,$5,"+"}' | gzip > ${tmp}/${newname}.bed.gz
   fi
   echo 'Number of mappable reads:' >> ${output}/${newname}.align.log
   readCount=`zcat ${tmp}/${newname}.bed.gz | grep "" -c`
   echo ${readCount} >> ${output}/${newname}.align.log
   ## Remove rRNA and reverse the strand (PRO-seq).
   zcat ${tmp}/${newname}.bed.gz | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > ${tmp}/${newname}.nr.rs.bed.gz
   echo 'Number of mappable reads (excluding rRNA):' >> ${output}/${newname}.align.log
   echo `zcat ${tmp}/${newname}.nr.rs.bed.gz | grep "" -c` >> ${output}/${newname}.align.log
   ## Convert to bedGraph ... Cannot gzip these, unfortunately.
   bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.bed.gz -g ${CHINFO} -strand + > ${tmp}/${newname}\_plus.bedGraph
   bedtools genomecov -bg -i ${tmp}/${newname}.nr.rs.bed.gz -g ${CHINFO} -strand - > ${tmp}/${newname}\_minus.bedGraph
   ## Invert minus strand.
#   cat ${tmp}/$newname\_minus.noinv.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,-1*$4}' > ${tmp}/$newname_minus.bedGraph ## Invert read counts on the minus strand.
   ## normalized by RPM
   cat ${tmp}/${newname}_plus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_plus.rpm.bedGraph
   cat ${tmp}/${newname}_minus.bedGraph | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4*1000*1000/'$readCount'/1}' > ${tmp}/${newname}_minus.rpm.bedGraph
   ## Then to bigWig (nomalized and non-nomrmalized ones)
   bedGraphToBigWig ${tmp}/${newname}_plus.rpm.bedGraph ${CHINFO} ${output}/${newname}_plus.rpm.bw
   bedGraphToBigWig ${tmp}/${newname}_minus.rpm.bedGraph ${CHINFO} ${output}/${newname}_minus.rpm.bw
   bedGraphToBigWig ${tmp}/${newname}_plus.bedGraph ${CHINFO} ${output}/${newname}_plus.bw
   bedGraphToBigWig ${tmp}/${newname}_minus.bedGraph ${CHINFO} ${output}/${newname}_minus.bw
 done
fi

echo "QC" > ${output}/proseq2.0_read_report_${tmp}.log
for name in ${files[*]}
  do
  cat ${tmp}/${name}.QC.log >> ${output}/proseq2.0_read_report_${tmp}.log
done
echo "" >> ${output}/proseq2.0_read_report_${tmp}.log
echo "Mapping" >> ${output}/proseq2.0_read_report_${tmp}.log
for f in ${tmp}/*.sort.bam
 do j=`echo $f | awk -F"/" '{print $NF}' | rev | cut -d \. -f 3- |rev `
 cat ${output}/${j}.align.log >> ${output}/proseq2.0_read_report_${tmp}.log
done
