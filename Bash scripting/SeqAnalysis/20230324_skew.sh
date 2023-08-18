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
#we'll skip g monomers, that really hasn't felt useful.
##################
## Calculations ##
##################
#echo "sorting bed"
#date
#sort -k1,1 -k2,2n minimal_set.bed | gzip > min_set_sort.bed.gz

#echo "flattening bed"
#date
#zcat min_set_sort.bed.gz | bedtools merge -i stdin | gzip > min_set_flat.bed.gz

#########
#echo "making windows"
#echo "making staggered wins"
#date
#zcat min_set_flat.bed.gz | bedops --chop 10 --stagger 1 -x - | gzip > stag.txt.gz

#echo "making disjoint wins"
#date
#zcat min_set_flat.bed.gz | bedops --chop 10 -x - | gzip > nstag.txt.gz

##########
#echo "Generating window content"
#date
#gen=/scratch/mford5/references/hg38/hg38.fa
#zcat stag.txt.gz | bedtools nuc -fi ${gen} -bed stdin -seq | gzip > stag_nuc.txt.gz
####################
#echo "removing headers"
#date
#zcat stag_nuc.txt.gz | sed '1d' | gzip > nh_stag_nuc.txt.gz 
####################
#echo "making GC content beds"
#echo "for staggered"
#date
#zcat nh_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, "ID", $5}' | gzip > stag_GC.txt.gz
####################
#echo "making GC skew beds"
#echo "for staggered"
#date
#zcat nh_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {
#	if ($5==0) {
#		print $1, $2, $3, "ID", "0"
#	} else {
#		dif=$8-$7; sum=$8+$7; skew=dif/sum; print $1, $2, $3, "ID", skew
#		}
#	}' | gzip > stag_skew.txt.gz
#####################
### aid apobec motif TCW (A/T) and WGA, -seq then check for any of those.###
echo "calc apobec"
date
zcat nh_stag_nuc.txt.gz | awk -F "\t" 'BEGIN {OFS=FS} {if ($13 ~ /^TC[AT]|^[AT]GA/) {
	print $1, $2, $2+3, "ID", "1"}}' | bedops --merge - | gzip > motifs.txt.gz
#basically this should make 3bp windows with a score of 1 if they start with the apobec motif, using the staggered windows I already have. I don't really want to make this into a bw, I don't think that's accurate representation.
##############################
echo "smoothing staggered GC"
date
gunzip nstag.txt.gz
zcat stag_GC.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean nstag.txt - > smooth_GC.txt
echo "smoothing staggered skew"
date
zcat stag_skew.txt.gz | bedmap --echo --delim "\t" --unmapped-val 0 --mean nstag.txt - > smooth_skew.txt
awk -F "\t" 'BEGIN {OFS=FS} {print $1, $2, $3, $4*-1}' smooth_skew.txt > m_smooth_skew.txt

echo "getting bigwigs"
date
bg2bw=/scratch/mford5/tools/bedGraphToBigWig
chrsizes=/scratch/mford5/references/hg38/hg38.chrom.sizes
${bg2bw} smooth_GC.txt ${chrsizes} smooth_GC.bw
${bg2bw} smooth_skew.txt ${chrsizes} p_smooth_skew.bw
${bg2bw} m_smooth_skew.txt ${chrsizes} m_smooth_skew.bw
echo "fin"
date
