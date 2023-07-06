#!/bin/bash
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=groseq_%j.out
#SBATCH --mail-type=all
#remember to change 1g back to 50 and 1 hr back to 50
module load bwa
module load cutadapt
module load samtools
module load bedtools
module load seqtk
module load fastx-toolkit
module load prinseq
module load bedops
module load perl/5.10.1/b1

bedGraphToBigWig () {
/scratch/mford5/tools/bedGraphToBigWig
}
export -f bedGraphToBigWig

#bwa index /scratch/mford5/references/hg38/hg38.fa
export bwaIndex=/scratch/mford5/references/hg38/hg38.fa
export chromInfo=/scratch/mford5/references/hg38/hg38.chrom.sizes
#bash /scratch/mford5/tools/proseq2.0-master/proseq2.0.bsh -i $bwaIndex -c $chromInfo -SE -G \
#-T tmp -O groseq_results -I MCF7.Kraus.Rep3.Groseq 

#pipeline fails on bedgraph to bigwig. IDK why - it's clear it calls the correct function, but it does not pass the input to bedgraph to bigwig correctly. As I have absolutely no idea
#why, I'm going to have to just brute force the final step through and spend some time trying to figure out wtf is going on later. Probably something to do with me redefining the function to get it to work in subshells - not sure if the new definition is correctly passing the input arguments into the real function.


/scratch/mford5/tools/bedGraphToBigWig MCF7.Kraus.Rep3.Groseq_QC_minus.rpm.bedGraph $chromInfo R3.groseq.minus.rpm.bw
/scratch/mford5/tools/bedGraphToBigWig MCF7.Kraus.Rep3.Groseq_QC_plus.rpm.bedGraph $chromInfo R3.groseq.plus.rpm.bw
/scratch/mford5/tools/bedGraphToBigWig MCF7.Kraus.Rep3.Groseq_QC_minus.bedGraph $chromInfo R3.groseq.minus.bw
/scratch/mford5/tools/bedGraphToBigWig MCF7.Kraus.Rep3.Groseq_QC_plus.bedGraph $chromInfo R3.groseq.plus.bw
