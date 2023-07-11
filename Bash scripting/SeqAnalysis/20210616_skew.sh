#!/bin/bash

#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output="import_%u_%j.out"
#SBATCH --mail-type=all

#./twoBitToFa hg38.2bit hg38.fa

#module load bedtools
#module load samtools
#module load deeptools

#bedtools makewindows -b region.sizesorted.bed -w 10 -s 5 > win_reg.txt #makes 10bp windows sliding 5 ea time across human genome.
#bedtools nuc -fi hg38.fa -bed win_reg.txt > totalcontent.txt #generates file for nucleotide content in each window. 
#awk -F "\t" 'BEGIN {OFS=FS} {if (($7+$8) != 0) {print}}' totalcontent.txt > cleaning.txt
#sed '1d' cleaning.txt > cleaning1.txt
#rm cleaning.txt
#awk -F "\t" '{print $0, ($7-$8), ($7+$8)}' OFS='\t' cleaning1.txt > cleaning2.txt
#rm cleaning1.txt
awk -F "\t" '{print $1, $2, $3, $13/$14}' OFS='\t' cleaning2.txt > skew.txt
#rm cleaning2.txt
#grep -v 'chrX' skew.txt > skew.tr1.txt
#grep -v 'chrY' skew.tr1.txt > skew.tr2.txt
#awk '{if (NR>296861242 && NR<296861246){print $0} }' skew.tr2.txt
#bedtools sort -i skew.tr2.bed
sort -k1,1 -k2,2n skew.txt >skew.sort.txt
./bedGraphToBigWig skew.sort.txt chromSizes.txt skew.bw
