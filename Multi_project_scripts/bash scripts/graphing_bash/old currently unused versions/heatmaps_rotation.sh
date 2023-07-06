#!/bin/bash
#SBATCH --mem=48G
#SBATCH --time=96:00:00
#SBATCH --output="troubleshoot/heatmap_pipeline_tst2.%j.out"
#SBATCH --mail-type=all

###Above tells SLURM this is a slurm script, that I want it to allocate 48G memory for, limit to running 96 hours, and email me when it starts/finishes/fails.

module load samtools
module load deeptools
for f in MCF7.nullFiles/*.bam;
       do
      	samtools index $f #generates bam.bai for every bam file in the folder.
	outFile=$(basename "$f" .bam | awk '{print "pipeOut/"$1".sorted.bam"}') #makes a file name for output. Takes filename - .bam, then adds .sorted.bam to that string.
        samtools sort -o $outFile #sorts the bam file by chromosome then order in chromosome.
        outFile2=$(basename "$f" .bam | awk '{print "pipeOut/"$1".bw"}') #makes filename by taking bam off filename then adding bw to it.
bamCoverage --bam MCF7.nullFiles/MCF7.1.null.sorted.bam \
        --outFileName MCF7.1.null.test2.sorted.bw \
        --binSize 10 \
	--normalizeUsingRPKM
done

#cd pipeOut

#bigwigCompare --bigwig1=MCF7.1.null.sorted.bw --bigwig2=MCF7.2.null.sorted.bw --outFileName=merge.tmp.bw --outFileFormat=bigwig
#bigwigCompare --bigwig1=merge.tmp.bw --bigwig2=MCF7.3.null.sorted.bw --outFileName=merge2.tmp.bw --outFileFormat=bigwig
#bigwigCompare --bigwig1=merge2.tmp.bw --bigwig2=MCF7.4.null.sorted.bw --outFileName=merge3.tmp.bw --outFileFormat=bigwig
#bigwigCompare --bigwig1=merge3.tmp.bw --bigwig2=MCF7.5.null.sorted.bw --outFileName=merge4.tmp.bw --outFileFormat=bigwig
#bigwigCompare --bigwig1=merge4.tmp.bw --bigwig2=MCF7.6.null.sorted.bw --outFileName=merge5.tmp.bw --outFileFormat=bigwig
#bigwigCompare --bigwig1=merge5.tmp.bw --bigwig2=MCF7.7.null.sorted.bw --outFileName=merge.bw --outFileFormat=bigwig
#there has to be an easier way to do this but unfortunately every iterative option I tried failed.
#rm *.temp.bw

computeMatrix reference-point -S=MCF7.1.null.test2.sorted.bw -R=cpg.bed -a=2000 -b=2000 --outFileName=test3.cpg.matrix --referencePoint=center
plotHeatmap -m=test3.cpg.matrix -out=kdm5.test3.heatmap.pdf --colorMap=coolwarm --refPointLabel=center #see deeptools issue 683. Version installe on bluehive ignores reference point label, instead labels center always as TSS, therefore manual labeling specified here.
