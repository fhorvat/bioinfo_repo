#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_to_genome.bismark
#PBS -l select=ncpus=20:mem=60g
#PBS -j oe
#PBS -J 0-%N_SAMPLES
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh

# override number of threads 
THREADS=4

# input
UNIQ_SEQ=($(printf "%s\n" "${IN_SEQ[@]%_[1,2,s].txt.gz}" | sort -u))
FILE=${UNIQ_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.txt.gz}

# mapping parameters
SCRIPT_PAR="--genome_folder ${BISMARK_INDEX} \
--parallel ${THREADS} \
--unmapped \
--output_dir . \
--temp_dir ."

# ----------------Commands------------------- #
# mapping
if [ ${SINGLE_END} == TRUE ]
then

   # run bismark
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} --single_end ${FILE}.txt.gz

   # rename
   mv ${BASE}.txt.gz_bismark_bt2.bam ${BASE}_bismark_bt2.bam
   mv ${BASE}.txt.gz_bismark_bt2_SE_report.txt ${BASE}_bismark_bt2_SE_report.txt
   mv ${BASE}.txt.gz_unmapped_reads.fq.gz ${BASE}.unmapped.txt.gz

   # sort bam
   samtools sort -@ ${THREADS} -o ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

   # rename over original file
   [ -f "${BASE}_bismark_bt2.sorted.bam" ] && mv ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2.bam

   # bam index
   samtools index -@ ${THREADS} ${BASE}_bismark_bt2.bam

   # bam to scaled bedGraph, bedGraph to bigWig
   genomeCoverageBed -ibam ${BASE}_bismark_bt2.bam -bg -split > ${BASE}.bedGraph
   wigToBigWig ${BASE}.bedGraph ${CHR_LENGTH} ${BASE}_bismark_bt2.raw.bw
   [ -f "${BASE}_bismark_bt2.raw.bw" ] && rm ${BASE}.bedGraph


else

   # run bismark
   bismark ${SCRIPT_PAR} --genome_folder ${BISMARK_INDEX} -1 ${FILE}_1.txt.gz -2 ${FILE}_2.txt.gz ${FILE}_s.txt.gz

   # rename
   mv ${BASE}_1.txt.gz_bismark_bt2_pe.bam ${BASE}_bismark_bt2_pe.bam
   mv ${BASE}_1.txt.gz_bismark_bt2_PE_report.txt ${BASE}_bismark_bt2_PE_report.txt
   mv ${BASE}_1.txt.gz_unmapped_reads_1.fq.gz ${BASE}_1.unmapped.txt.gz
   mv ${BASE}_2.txt.gz_unmapped_reads_2.fq.gz ${BASE}_2.unmapped.txt.gz

   # sort bam
   samtools sort -@ ${THREADS} -o ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2_pe.bam

   # rename over original file
   [ -f "${BASE}_bismark_bt2.sorted.bam" ] && mv ${BASE}_bismark_bt2.sorted.bam ${BASE}_bismark_bt2_pe.bam

   # bam index
   samtools index -@ ${THREADS} ${BASE}_bismark_bt2_pe.bam

   # bam to scaled bedGraph, bedGraph to bigWig
   genomeCoverageBed -ibam ${BASE}_bismark_bt2_pe.bam -bg -split > ${BASE}.bedGraph
   wigToBigWig ${BASE}.bedGraph ${CHR_LENGTH} ${BASE}_bismark_bt2_pe.raw.bw
   [ -f "${BASE}_bismark_bt2_pe.raw.bw" ] && rm ${BASE}.bedGraph

fi
