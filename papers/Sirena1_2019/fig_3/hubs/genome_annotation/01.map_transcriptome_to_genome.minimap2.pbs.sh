#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=50g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.map_transcriptome_to_genome.minimap2
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=50G

IN_DIR=..
IN_SEQ=(`find ${IN_DIR} -name "*.fa"`)
FILE=${IN_SEQ[0]}
BASE=${FILE#${IN_DIR}/}
BASE=${BASE%.fa}

FILE_TRANS=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20190702.cdna.all.fa.gz

MINIMAP2_PAR="-x splice \
-a \
-t $THREADS"

# ----------------Commands------------------- #
# align assembled transcriptome to genome
minimap2 $MINIMAP2_PAR $FILE $FILE_TRANS | \
sambamba view -f bam -F "not (unmapped or supplementary)" -t 6 -S /dev/stdin | \
sambamba sort -t 6 -o ${BASE}.mm10.ensembl.93.cdna.bam /dev/stdin
