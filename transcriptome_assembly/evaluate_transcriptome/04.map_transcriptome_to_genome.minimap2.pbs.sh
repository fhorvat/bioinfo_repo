#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.04.map_transcriptome_to_genome.minimap2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

BASE_DIR=`basename ${PWD/trinity*/}`

IN_DIR_GEN=../../../genomes/$BASE_DIR
IN_SEQ_GEN=(${IN_DIR_GEN}/*.fa.gz)
FILE_GEN=${IN_SEQ_GEN[0]}

IN_DIR_TRANS=.
IN_SEQ_TRANS=(${IN_DIR_TRANS}/*.fasta)
FILE_TRANS=${IN_SEQ_TRANS[0]}

BASE=${FILE_TRANS#${IN_DIR_TRANS}/}
BASE=${BASE/.fasta/.t2g}

MINIMAP2_PAR="-x splice \
-a \
-t $THREADS"

# ----------------Commands------------------- #
# align assembled transcriptome to genome
minimap2 $MINIMAP2_PAR $FILE_GEN $FILE_TRANS | \
sambamba view -f bam -F "not (unmapped or supplementary)" -t 6 -S /dev/stdin | \
sambamba sort -t 6 -o ${BASE}.bam /dev/stdin

# count mapped reads
sambamba view -t $THREADS ${BASE}.bam | awk '{print $1}' | sort | uniq | wc -l > ${BASE}.read_counts.txt

# count all reads in fasta file
grep -c ">" $FILE_TRANS >> ${BASE}.read_counts.txt
