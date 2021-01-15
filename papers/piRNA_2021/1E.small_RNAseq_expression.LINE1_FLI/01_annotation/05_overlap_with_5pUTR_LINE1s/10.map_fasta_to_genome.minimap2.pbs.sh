#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=20g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.map_fasta_to_genome.minimap2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=20G

IN_DIR_GEN=/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed
IN_SEQ_GEN=($(find ${IN_DIR_GEN} -maxdepth 1 -name "*.fasta"))
FILE_GEN=${IN_SEQ_GEN[0]}

IN_DIR_TRANS=.
IN_SEQ_TRANS=($(find ${IN_DIR_TRANS} -maxdepth 1 -name "*.no_gaps.fasta"))
FILE_TRANS=${IN_SEQ_TRANS[0]}

BASE=${FILE_TRANS#${IN_DIR_TRANS}/}
BASE=${BASE%.fasta}

MINIMAP2_PAR="-x asm5 \
-a \
-t $THREADS"

# ----------------Commands------------------- #
# align assembled transcriptome to genome
minimap2 $MINIMAP2_PAR $FILE_GEN $FILE_TRANS | \
sambamba view -f bam -F "not (unmapped or supplementary)" -t $THREADS -S /dev/stdin | \
sambamba sort -t $THREADS -o ${BASE}.bam /dev/stdin
