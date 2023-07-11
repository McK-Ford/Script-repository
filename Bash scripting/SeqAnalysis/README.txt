2021 06 16 is bedtools strategy, across entire genome relative to the plus strand, using 5bp windows.

2021 07 19 is sense skew using bedtools strategy. Bad strategy really as it's meant to be sense by gene so you can't look at the bw and know the orientation of a given region and it doesn't work for overlapping regions.

2021 09 13 is skew, using BEDOPs strategy, with minus and plus calculated separately.

2023 03 17 is skew, GC percentage, and G monomer stretches based on refseq curated, using BEDOPs strategy, with minus and plus calculated separately. This is the script for the 20230317 data in the projects folder.

2023 03 24 skew is skew and GC percentage for extended windows (5kb around gencode all genes/ islands), using BEDOPs strategy, and with minus strand just calculated by inverting plus strand score. This is the script for the 20230329 data in the projects folder.