#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
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
echo "this file is used to generate skew in the sense direction for"
echo "a mixed plus minus region set. Be aware there will be uncovered regions"
echo "and that this skew bigwig is not generically applicable"
echo "Merging any overlapping regions in input beds to one big region"
bedops --merge gene_1k_win_p.bed > clean_win_p.bed
bedops --merge gene_1k_win_m.bed > clean_win_m.bed
echo "Making 5 bp windows across all bed regions" 
bedtools makewindows -b clean_win_p.bed -w 5 | gzip -f >pwin.txt.gz
bedtools makewindows -b clean_win_m.bed -w 5 | gzip -f >mwin.txt.gz

echo "Generating window bp content"
zcat pwin.txt.gz | bedtools nuc -fi /scratch/mford5/references/hg38/hg38.fa -bed stdin | \
	gzip -f > totalcontentp.txt.gz #generates file for nucleotide content in each window.
zcat mwin.txt.gz | bedtools nuc -fi /scratch/mford5/references/hg38/hg38.fa -bed stdin | \
	gzip -f > totalcontentm.txt.gz  

echo "removing headers for computational ease"
zcat totalcontentp.txt.gz | sed '1d' | gzip -f > no_headerp.txt.gz
zcat totalcontentm.txt.gz | sed '1d' | gzip -f > no_headerm.txt.gz

echo "doing the skew math"
echo "skew here is defined as (G-C)/(G+C)"
zcat no_headerp.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if (($8+$7)==0) {
	print $1, $2, $3, "0"} else {
	dif=$8-$7; sum=$8+$7; p_skew=dif/sum; print $1, $2, $3, p_skew}}' | \
	gzip -f > skewp.txt.gz
zcat no_headerm.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if (($8+$7) ==0) {
	print $1, $2, $3, "0"} else {
	dif=$7-$8; sum=$8+$7; m_skew=dif/sum; print $1, $2, $3, m_skew}}' | \
	gzip -f >skewm.txt.gz

echo "fixing the skew sort"
zcat skewp.txt.gz | sort-bed - --max-mem 12G > skewp.sort.txt
zcat skewm.txt.gz | sort-bed - --max-mem 12G > skewm.sort.txt

echo "separating out overlapping and nonoverlapping elements."
bedops --element-of skewp.sort.txt skewm.sort.txt > overlap1.bed
bedops --element-of skewm.sort.txt skewp.sort.txt > overlap2.bed
bedops --not-element-of 1 skewm.sort.txt skewp.sort.txt > nooverlap1.bed
bedops --not-element-of 1 skewp.sort.txt skewm.sort.txt > nooverlap2.bed
bedops --everything nooverlap1.bed nooverlap2.bed > merged.bed
echo "making the sense skew bigwig."
/scratch/mford5/tools/bedGraphToBigWig merged.bed /scratch/mford5/references/hg38/hg38.chrom.sizes senseskew_small.bw
