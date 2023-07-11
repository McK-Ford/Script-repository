peakstuff
checkingDeSeqBasedNormalization showed no change in VER0046 between kdm5 kd and no kd.
20220309_DiffBind_H3K27me3.R was because when we looked in IGV we seemed to notice some changes in H3K27me3 after KDM5B KO. I think I only identified a handful of genes that were significantly different, though.

TSS_problem
tsscallingnetcage.R :: Intersected NET-CAGE MCF7 data (GSE118075) with CpG islands, then attached these to the closest gencode TSSes to create a verified TSS dataset.

20211110 TruTSS stuff was an attempt to define the TSS by paired end pro-seq. Didn't work, 5' end never really piled up anywhere other than 50 bp before first pause peak.

shortest tss is most downstream that still intersects with CGI. Still biased but in opposite of all the longest TSS stuff I've done. Not really an adequate solution to the TSS problem.

annotation
seqTools.annotBedtoTxDbGene_ALTERED.r :: simple annotation function based off ben's and josh's old stuff.

PI_related
testing_BRGenomics.R :: Danko recommended BRGenomics over GroHMM for looking at pausing, bc GroHMM is not compatible with paired-end data. It didn't do exactly what I needed to do, but it provided a lot of valuable insights for the development of my own metaplotter.

Stranding_CGIs_by_transcription_05.19.2022.R :: Used GSE93229 to assign directions to CGI based on whether plus or minus had more tx. Goal was to make an annotation-agnostic way of looking at CGIs that still considered primary direction of tx. Also, would prove useful for looking at tx behaviour at enhancers in case any had preference. Also took location of maximal tx to act as a 'TSS'.

Strand_imputing_CGIs_06.07.2022.R :: Different way of assigning strands to CGI with a few more summary statistics than the 05.19.2022 attempt.

20220419_Classing_CGI_TSSes_By_expression was an attempt to bin the GSE93229 data by expression, then heatmap it. 

20220404_Calling_Pausing_Classes_new_fn is an attempt to call prox dist and silent based on the Danko data.

20220329_Refseq_Curated_CGI is just a refseq/cgi intersection filtering.
20220106_classing_danko_proseq.R classing genes by the danko data again
20220105_groseq_index.R I don't know why I called this GRO-seq because it's clearly using PRO-seq to class, though I did use shortest TSS here not longest.

TESTING, OLD, OTHERWISE NOT TO BE REUSED/ Graphing related

20220613_testingproseq.R :: GSE93229 proseq metaplotting attempt. Was trying to rework my metaplotter for using the stranded bams - I believe I ran into problems because the metaplotter script at the time expected reads on both strands.

20210607_color_CGIlen_heatmaps.R :: VER00046 cut&tag heatmaps, using an old version of my metaplotter. At promoter CGIs vs non-CGI promoters.

20220504FragnInset_testset1 :: I took the CGI TSSes, and metaplotted the 1 kb around the TSS for H3K4me3, H2AZ, and KDM5B with and without insert, with a filter applied to insert size, and with it resized to just the center 50 reads.

20220607Heatmaps based on strand_imputed_CGIs :: pretty self-explanatory. Uses the strand imputed CGIs from Strand_imputing_CGIs_06.07.2022.R.

20220520BnW_CGI_len_sorted_heatmaps :: refseq curated TSS VER00046 heatmaps.
20220511_making_cgi_len_heatmaps.R :: more refseq curated TSS VER00046 heatmaps.
20220510for_zach.R :: presentation-quality graphs to give to zach of refseq curated TSS VER00046 heatmaps and metaplots for a presentation.
20220413Intersecting_n_not_CGI_metaplots.r :: refseq curated longest TSS metaplots, metaplots of VER00046 with it split based on whether the genes intersected CpG islands or not.
20220414heatmapping.R :: more VER00046 heatmaps using a slightly different method.
20220414Scaled_Metaplots_new.R :: again, more VER00046 with refseq curated longest, though these ones are metaplots not heatmaps.
20220504 FragnInsert testsets 2 and 3 were I believe metaplots testing how filtering for different fragment lengths and whether or not to include the insert affected cut and tag plots.
Clean graphs for chromatin club is exactly as it sounds. Scaled metaplots new2 is another metaplotter, because apparently I decided when I first started that I needed to make a new graphing script every time I graphed anything.
also metaplotting_genes_example.R
