I made a lot of 'graphing' scripts that just used different files or slightly different variables and aren't really worth uploading here. I've kept the most complicated ones walking through special file manipulation, but here are the basics. Ultimately, I settled on a command line interaction for small scale graphing. (ought to write a script for a larger one, but some things like setting the z score just need to be done manually).

Can list files ie like so matrices=(mat1 mat2 mat3) and do for m in ${matrices[*]}


Format: 
module load deeptools
computeMatrix reference-point -S sample.bw sample2.bw -R regions.bed regions2.bed -a=downDist -b=upDist --outFileName=sample.matrix.gz --referencePoint=TSS
#common variant --referencePoint=center

plotHeatmap -m sample.matrix -out=sample.pdf --colorMap=someColorName --sortRegions=descend --sortUsing=region_length
#common variant sortRegions=keep
#common addendum --kmeans=3
--labelRotation=45
--outFileSortedRegions=tabs/${file}.reg.txt

#For scale, scale-regions instead of reference point, also need -m and -bs. Deeptools scaleregions scales everything to number of bins instead of bp.

To do this stranded, use a stranded bed (bed6), sorted by strand, and run it for both the + and - strand bigwigs.
Next, zcat sample.matrix.gz | sed "1d" > noheader.sample.txt for both plus and minus.
Then, get a header too. Need to modify header manually to account for new number in each section.
Use sed to split the parental file by sense and antisense data. For an example, the genes classed by the 20230329 skew, we use: sed -n "1,3606p" for prox plus, "3607,7073p" for prox minus, "7074,10301p" for dist plus, and "10302,13542p" for dist minus.
cat the plus and minus together.
To sort by skew or other external info instead of something within the mat like region size:
awk 'FNR==NR{id[$4] = $0;next};{$1 in id}{print id[$1]}' tmp/prox.tmp bed/prox_ord.bed > orderedmat.tmp
Cat the header back on and gzip
Then you can plot heatmap.