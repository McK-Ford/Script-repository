#!/bin/bash
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output="mk_skew_%j.out"

#####################
## Loading modules ##
#####################
module load bedtools
module load bedops

##################
## Calculations ##
##################
echo "sorting bed"
date
sort-bed plus_win_for_skew.bed --max-mem 64G | gzip > p_sort.bed.gz
sort-bed minus_win_for_skew.bed --max-mem 64G | gzip > m_sort.bed.gz

#########
echo "making windows"
echo "making staggered wins"
date
zcat p_sort.bed.gz | bedops --chop 10 --stagger 1 -x - | gzip > pwin_stag.txt.gz
zcat m_sort.bed.gz | bedops --chop 10 --stagger 1 -x - | gzip > mwin_stag.txt.gz

echo "making disjoint wins"
date
zcat p_sort.bed.gz | bedops --chop 10 -x - | gzip > pwin_nstag.txt.gz
zcat m_sort.bed.gz | bedops --chop 10 -x - | gzip > mwin_nstag.txt.gz

echo "making single bp wins (for G monomer len)"
date
zcat p_sort.bed.gz | bedops --chop 1 -x - | gzip > pwin_1bp.txt.gz
zcat m_sort.bed.gz | bedops --chop 1 -x - | gzip > mwin_1bp.txt.gz

##########
echo "Generating window content"
echo "for staggered wins"
date
gen=/scratch/mford5/references/hg38/hg38.fa
zcat pwin_stag.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > p_stag_nuc.txt.gz
zcat mwin_stag.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > m_stag_nuc.txt.gz

echo "for disjoint wins"
date
zcat pwin_nstag.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > p_nstag_nuc.txt.gz
zcat mwin_nstag.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > m_nstag_nuc.txt.gz

echo "for sbp wins (G monomer len)"
date
zcat pwin_1bp.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > p_1bp_nuc.txt.gz
zcat mwin_1bp.txt.gz | bedtools nuc -fi ${gen} -bed stdin | gzip > m_1bp_nuc.txt.gz

####################
echo "removing headers"
date
zcat p_stag_nuc.txt.gz | sed '1d' | gzip > nhp_stag_nuc.txt.gz 
zcat m_stag_nuc.txt.gz | sed '1d' | gzip > nhm_stag_nuc.txt.gz
zcat p_nstag_nuc.txt.gz | sed '1d' | gzip > nhp_nstag_nuc.txt.gz
zcat m_nstag_nuc.txt.gz | sed '1d' | gzip > nhm_nstag_nuc.txt.gz
zcat p_1bp_nuc.txt.gz | sed '1d' | gzip > nhp_1bp_nuc.txt.gz
zcat m_1bp_nuc.txt.gz | sed '1d' | gzip > nhm_1bp_nuc.txt.gz
####################
echo "making GC content beds"
echo "for staggered"
date
use ID as name placeholder
zcat nhp_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, "ID", $5}' | gzip > p_stag_GC.txt.gz
zcat nhm_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, "ID", $5}' | gzip > m_stag_GC.txt.gz
echo "for disjoint"
date
zcat nhp_nstag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, $5}' > p_nstag_GC.txt
zcat nhm_nstag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, $5}' > m_nstag_GC.txt
####################
echo "making GC skew beds"
echo "for staggered"
date
zcat nhp_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($5==0) {
	print $1, $2, $3, "ID", "0"} else {
	dif=$8-$7; sum=$8+$7; skew=dif/sum; print $1, $2, $3, "ID", skew}}' | gzip > p_stag_skew.txt.gz
zcat nhm_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($5==0) {
	print $1, $2, $3, "ID", "0"} else {
	dif=$7-$8; sum=$8+$7; skew=dif/sum; print $1, $2, $3, "ID", skew}}' | gzip > m_stag_skew.txt.gz
echo "for disjoint"
date
zcat nhp_nstag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($5==0) {
	print $1, $2, $3, "0"} else {
	dif=$8-$7; sum=$8+$7; skew=dif/sum; print $1, $2, $3, skew}}' > p_nstag_skew.txt
zcat nhm_nstag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($5==0) {
	print $1, $2, $3, "0"} else {
	dif=$7-$8; sum=$8+$7; skew=dif/sum; print $1, $2, $3, skew}}' > m_nstag_skew.txt #C skew on plus strand = G skew on minus strand
#####################
echo "calc monomer lens"
echo "getting only Gs and collapsing them"
date
zcat nhp_1bp_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($8>0) {
	print $1, $2, $3}}' | bedops --merge - | gzip > p_Gs.txt.gz
zcat nhm_1bp_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($7>0) {
	print $1, $2, $3}}' | bedops --merge - | gzip > m_Gs.txt.gz #C on plus strand = G on minus strand
echo "calc win len"
date
zcat p_Gs.txt.gz | awk -F "\t" '{if (($3-$2)>1){rlen=$3-$2-1; print $1, $2, $3, "ID", rlen}}' | gzip > p_Gwins.txt.gz #no credit for single bp Gs
zcat m_Gs.txt.gz | awk -F "\t" '{if (($3-$2)>1){rlen=$3-$2-1; print $1, $2, $3, "ID", rlen}}' | gzip > m_Gwins.txt.gz
##############################
echo "smoothing staggered GC"
date
gunzip pwin_nstag.txt.gz
gunzip mwin_nstag.txt.gz
zcat p_stag_GC.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean pwin_nstag.txt - > p_smooth_GC.txt
zcat m_stag_GC.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean mwin_nstag.txt - > m_smooth_GC.txt
echo "smoothing staggered skew"
date
zcat p_stag_skew.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean pwin_nstag.txt - > p_smooth_skew.txt
zcat m_stag_skew.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean mwin_nstag.txt - > m_smooth_skew.txt
echo "getting windowed monomer graph"
echo "mean"
date
zcat p_Gwins.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean pwin_nstag.txt - > p_Gmono_mean.txt
#zcat m_Gwins.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean mwin_nstag.txt - > m_Gmono_mean.txt
echo "max"
zcat p_Gwins.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --max pwin_nstag.txt - > p_Gmono_max.txt
#zcat m_Gwins.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --max mwin_nstag.txt - > m_Gmono_max.txt

echo "getting bigwigs"
date
bg2bw=/scratch/mford5/tools/bedGraphToBigWig
chrsizes=/scratch/mford5/references/hg38/hg38.chrom.sizes
${bg2bw} p_nstag_GC.txt ${chrsizes} p_nstag_GC.bw
${bg2bw} m_nstag_GC.txt ${chrsizes} m_nstag_GC.bw
${bg2bw} p_nstag_skew.txt ${chrsizes} p_nstag_skew.bw
${bg2bw} m_nstag_skew.txt ${chrsizes} m_nstag_skew.bw
${bg2bw} p_smooth_GC.txt ${chrsizes} p_smooth_GC.bw
${bg2bw} m_smooth_GC.txt ${chrsizes} m_smooth_GC.bw
${bg2bw} p_smooth_skew.txt ${chrsizes} p_smooth_skew.bw
${bg2bw} m_smooth_skew.txt ${chrsizes} m_smooth_skew.bw
${bg2bw} p_Gmono_mean.txt ${chrsizes} p_Gmono_mean.bw
${bg2bw} m_Gmono_mean.txt ${chrsizes} m_Gmono_mean.bw
${bg2bw} p_Gmono_max.txt ${chrsizes} p_Gmono_max.bw
${bg2bw} m_Gmono_max.txt ${chrsizes} m_Gmono_max.bw
echo "fin"
date

#echo "fixing the skew sort"
#zcat skewp.txt.gz | sort-bed - --max-mem 30G > skewp.sort.txt
#zcat skewm.txt.gz | sort-bed - --max-mem 30G > skewm.sort.txt
#echo "making directional bigwigs"
#/scratch/mford5/tools/bedGraphToBigWig skewp.sort.txt \
#	/scratch/mford5/references/hg38/hg38.chrom.sizes plus_strand_skew_10.bw
#/scratch/mford5/tools/bedGraphToBigWig skewm.sort.txt \
#	/scratch/mford5/references/hg38/hg38.chrom.sizes minus_strand_skew_10.bw
#output is staggered (smoothed) GC percent and skew, disjoint GC percent and skew, and max and average monomer score for bins across each reference element.
