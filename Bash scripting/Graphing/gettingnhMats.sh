#!/bin/bash
#SBATCH -p standard -t 8:00:00
#SBATCH --output="matsnh%j.out"
#SBATCH -c 4 --mem=16G

mkdir noheader
mkdir tmp
mkdir processed_mats
mkdir pdf

module load deeptools
declare -a not_stranded=("apobec_motifs" "GC" "H3K79me2" "S9.6" "ttseq2")
declare -a stranded=("proseq" "skew" "ttseq")

for FILE in genmats/*.matrix.gz
do
prefix=$(basename ${FILE} .matrix.gz)
zcat genmats/${prefix}.matrix.gz | sed "1d" > noheader/${prefix}.noheader.txt
done

for FILE in tss_mats/*.matrix.gz
do
prefix=$(basename ${FILE} .matrix.gz)
zcat tss_mats/${prefix}.matrix.gz | sed "1d" > noheader/${prefix}.noheader.txt
done

for FILE in cgi_mats/*.matrix.gz
do
prefix=$(basename ${FILE} .matrix.gz)
zcat cgi_mats/${prefix}.matrix.gz | sed "1d" > noheader/${prefix}.noheader.txt
done

for FILE in "${not_stranded[@]}"
do
#no splicing or ordering treatment needed for basic file.
cat MODDED.header.txt noheader/${FILE}.noheader.txt | gzip > processed_mats/${FILE}.mat.gz
plotHeatmap -m processed_mats/${FILE}.mat.gz -o pdf/${FILE}.pdf --sortUsing=region_length

sed -n "1,3606p" noheader/${FILE}_tss.noheader.txt > tmp/prox_plus.tmp
sed -n "3606,7073p" noheader/${FILE}_tss.noheader.txt > tmp/prox_minus.tmp
sed -n "7073,10301p" noheader/${FILE}_tss.noheader.txt > tmp/dist_plus.tmp
sed -n "10302,13542p" noheader/${FILE}_tss.noheader.txt > tmp/dist_minus.tmp
cat tmp/prox_plus.tmp tmp/prox_minus.tmp > tmp/prox.tmp
cat tmp/dist_plus.tmp tmp/dist_minus.tmp > tmp/dist.tmp

#prox order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_tss_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_tss_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_tss_proxskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_tss_proxskew.mat.gz -o pdf/${FILE}_tss_proxskew.pdf --sortRegions=no

#dist order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_cgi_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_cgi_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_tss_distskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_tss_distskew.mat.gz -o pdf/${FILE}_tss_distskew.pdf --sortRegions=no

sed -n "1,3606p" noheader/${FILE}_cgi.noheader.txt > tmp/prox_plus.tmp
sed -n "3607,7073p" noheader/${FILE}_cgi.noheader.txt > tmp/prox_minus.tmp
sed -n "7074,10301p" noheader/${FILE}_cgi.noheader.txt > tmp/dist_plus.tmp
sed -n "10302,13542p" noheader/${FILE}_cgi.noheader.txt > tmp/dist_minus.tmp
cat tmp/prox_plus.tmp tmp/prox_minus.tmp > tmp/prox.tmp
cat tmp/dist_plus.tmp tmp/dist_minus.tmp > tmp/dist.tmp
#prox order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_tss_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_tss_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_cgi_proxskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_cgi_proxskew.mat.gz -o pdf/${FILE}_cgi_proxskew.pdf --sortRegions=no

#dist order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_cgi_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_cgi_order.bed > tmp/dist_ord.tmp
cat header_cgi_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_cgi_distskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_cgi_distskew.mat.gz -o pdf/${FILE}_cgi_distskew.pdf --sortRegions=no
done

for FILE in "${stranded[@]}"
do

sed -n "1,3606p" noheader/${FILE}_plus.noheader.txt > tmp/prox_plus.tmp
sed -n "1,3467p" noheader/${FILE}_minus.noheader.txt > tmp/prox_minus.tmp
sed -n "3607,6834p" noheader/${FILE}_plus.noheader.txt > tmp/dist_plus.tmp
sed -n "3468,6708p" noheader/${FILE}_minus.noheader.txt > tmp/dist_minus.tmp
#no ordering needed for basic file
cat MODDED.header.txt tmp/prox_plus.tmp tmp/prox_minus.tmp tmp/dist_plus.tmp tmp/dist_minus.tmp | gzip > processed_mats/${FILE}.mat.gz
plotHeatmap -m processed_mats/${FILE}.mat.gz -o pdf/${FILE}.pdf --sortUsing=region_length

sed -n "1,3606p" noheader/${FILE}_tss_plus.noheader.txt > tmp/prox_plus.tmp
sed -n "1,3467p" noheader/${FILE}_tss_minus.noheader.txt > tmp/prox_minus.tmp
sed -n "3607,6834p" noheader/${FILE}_tss_plus.noheader.txt > tmp/dist_plus.tmp
sed -n "3468,6708p" noheader/${FILE}_tss_minus.noheader.txt > tmp/dist_minus.tmp
cat tmp/prox_plus.tmp tmp/prox_minus.tmp > tmp/prox.tmp
cat tmp/dist_plus.tmp tmp/dist_minus.tmp > tmp/dist.tmp


#prox order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_tss_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_tss_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_tss_proxskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_tss_proxskew.mat.gz -o pdf/${FILE}_tss_proxskew.pdf --sortRegions=no

#dist order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_cgi_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_cgi_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_tss_distskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_tss_distskew.mat.gz -o pdf/${FILE}_tss_distskew.pdf --sortRegions=no

sed -n "1,3606p" noheader/${FILE}_cgi_plus.noheader.txt > tmp/prox_plus.tmp
sed -n "1,3467p" noheader/${FILE}_cgi_minus.noheader.txt > tmp/prox_minus.tmp
sed -n "3607,6834p" noheader/${FILE}_cgi_plus.noheader.txt > tmp/dist_plus.tmp
sed -n "3468,6708p" noheader/${FILE}_cgi_minus.noheader.txt > tmp/dist_minus.tmp
cat tmp/prox_plus.tmp tmp/prox_minus.tmp > tmp/prox.tmp
cat tmp/dist_plus.tmp tmp/dist_minus.tmp > tmp/dist.tmp

#prox order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_tss_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_tss_order.bed > tmp/dist_ord.tmp
cat header_tss_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_cgi_proxskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_cgi_proxskew.mat.gz -o pdf/${FILE}_cgi_proxskew.pdf --sortRegions=no

#dist order
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_cgi_order.bed > tmp/prox_ord.tmp
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/dist.tmp bed/dist_cgi_order.bed > tmp/dist_ord.tmp
cat header_cgi_modded.txt tmp/prox_ord.tmp tmp/dist_ord.tmp | gzip > processed_mats/${FILE}_cgi_distskew.mat.gz
plotHeatmap -m processed_mats/${FILE}_cgi_distskew.mat.gz -o pdf/${FILE}_cgi_distskew.pdf --sortRegions=no

done
