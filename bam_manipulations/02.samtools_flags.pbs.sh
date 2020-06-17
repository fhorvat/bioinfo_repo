#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N samtools_count
#PBS -l select=ncpus=1:mem=10g
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #

# ----------------Commands------------------- #
for bam_file in `ls *.bam`
do
echo -e ${bam_file%_multimap.bam} "\t" `samtools view -F 0x4 $bam_file | cut -f 1 | sort | uniq | wc -l` >> multimap_stats.log
done
