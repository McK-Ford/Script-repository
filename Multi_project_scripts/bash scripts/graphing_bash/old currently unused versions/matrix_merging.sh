#!/bin/bash
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=matrices_%j.out
#SBATCH --mail-type=all

mat=(R3_rpm_groseq groseq_rep3 skew_10)
for i in {0..2};
	do 
	zcat matrices/${mat[${i}]}_plus.matrix | awk 'NR<=2998 {print $0}' > Dist_p.tmp
	zcat matrices/${mat[${i}]}_plus.matrix | awk 'NR>=2998 {print $0}' | sed "1d" > TSS_p.tmp
	zcat matrices/${mat[${i}]}_minus.matrix | awk 'NR<=2815 {print $0}' | sed "1d" > Dist_m.tmp
	zcat matrices/${mat[${i}]}_minus.matrix | awk 'NR>=2815 {print $0}' | sed "1d" > TSS_m.tmp
	cat Dist_p.tmp Dist_m.tmp > Dist.tmp
	cat TSS_p.tmp TSS_m.tmp > TSS.tmp
	cat Dist.tmp TSS.tmp | gzip > ${mat[${i}]}.matrix
	rm *.tmp
done
