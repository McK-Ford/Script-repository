This went through a lot of versions, whether it was due to meaningful differences or due to lots of troubleshooting. Tried to distill it down to the most important ones.

2022 11 22 Julius Bowtie is the heavily modified version of the Julius pipeline I currently use, and all future edits should be off this template.

2021 11 30 much earlier version of julius pipeline. I no longer use it.

2021 11 11 Danko BWA I believe was created due to me realizing some of the flaws in 10 11 script, specifically for PE as that's what the edits pertain to.

2021 10 11 Danko BWA - Danko method, using BWA. Our version was designed for both PE and SE, but I believe it had some flaws. Not refsheet based unlike my later work. Processed by cutadapt instead of fastp. Uses perl prinseq for UMI management. Minimal QC compared to current version. Additionally, am unsure if this suffers from the paired end problem we noted with the julius pipeline where it wasn't built for PE data so it mapped PE reads to both strands (as second in pair is opposite strand). Ultimately did not use this further.