#!/bin/bash
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output="mk_skew_%j.out"
#SBATCH --mail-type=all

#####################
## Loading modules ##
#####################
module load bedtools
module load bedops

########################
## Making Skew bigwig ##
########################
sort-bed plus_win_for_skew.bed --max-mem 30G > p.sort.bed
sort-bed minus_win_for_skew.bed --max-mem 30G > m.sort.bed

bedops --chop 10 -x p.sort.bed | gzip >pwin.txt.gz
bedops --chop 10 -x m.sort.bed | gzip >mwin.txt.gz

#echo "Generating window content"
zcat pwin.txt.gz | bedtools nuc -fi /scratch/mford5/references/hg38/hg38.fa -bed stdin | \
	gzip > totalcontentp.txt.gz #generates file for nucleotide content in each window.
zcat mwin.txt.gz | bedtools nuc -fi /scratch/mford5/references/hg38/hg38.fa -bed stdin | \
        gzip > totalcontentm.txt.gz

#echo "removing headers"
zcat totalcontentp.txt.gz | sed '1d' | gzip > no_headerp.txt.gz
zcat totalcontentm.txt.gz | sed '1d' | gzip > no_headerm.txt.gz

#echo "doing the skew math"
zcat no_headerp.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if (($8+$7)==0) {
	print $1, $2, $3, "0"} else {
	dif=$8-$7; sum=$8+$7; p_skew=dif/sum; print $1, $2, $3, p_skew}}' | \
	gzip >skewp.txt.gz
zcat no_headerm.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if (($8+$7) ==0) {
	print $1, $2, $3, "0"} else {
	dif=$7-$8; sum=$8+$7; m_skew=dif/sum; print $1, $2, $3, m_skew}}' | \
	gzip >skewm.txt.gz

echo "fixing the skew sort"
zcat skewp.txt.gz | sort-bed - --max-mem 30G > skewp.sort.txt
zcat skewm.txt.gz | sort-bed - --max-mem 30G > skewm.sort.txt

#bedops --merge skewp.sort.txt > skewp.collapsed.txt
#bedops --merge skewm.sort.txt > skewm.collapsed.txt

echo "making directional bigwigs"
/scratch/mford5/tools/bedGraphToBigWig skewp.sort.txt \
	/scratch/mford5/references/hg38/hg38.chrom.sizes plus_strand_skew_10.bw
/scratch/mford5/tools/bedGraphToBigWig skewm.sort.txt \
	/scratch/mford5/references/hg38/hg38.chrom.sizes minus_strand_skew_10.bw

