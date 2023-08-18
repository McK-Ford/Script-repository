#load in libs
library(BRGenomics) # local install
library(GenomicRanges) #circ installed
library(Rsamtools) #circ installed
library(GenomicAlignments) #circ installed

#get cores
tmpvar <- Sys.getenv("SLURM_CPUS_PER_TASK")
ncores <- if (tmpvar == "") 1 else as.integer(tmpvar)
print(paste0("number of cores is: ", ncores))

####read bam files into GRanges list objects
sampstmp <- read.table("samples.txt")
samps <- as.list(sampstmp[[1]])
print(paste0("samples are ", samps))

getBams <- function(tmpf){
	bam_gr = GRanges(readGAlignmentPairs(paste0("bam/", tmpf, ".sorted.dupMark.bam")))
	print(head(bam_gr))
	print(tail(bam_gr))
	print(seqnames(bam_gr))
	return(bam_gr)
}

tmpL <- lapply(samps, getBams)

#no need to reinvent the wheel if this can do it
SampGRL <- as(tmpL, "GRangesList")
count_tab <- getSpikeInCounts(SampGRL, si_pattern = "spike", ncores = ncores, field = NULL)
count_tab$spikeRPM <- count_tab$spike_reads/1000000
write.table(count_tab, "full_count_tab.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(count_tab$spikeRPM, "factors_tab.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
