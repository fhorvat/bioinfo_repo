#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=12:mem=60g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.annotate_genome.Liftoff
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=12
MEMORY=60G

# target genome
INPUT_DIR=/common/DB/genome_reference/Muridae/Acomys_cahirinus/AcoCah.GCA_004027535.1
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "GCA_004027535.1_AcoCah_v1_BIUU.fa"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fa}

# query genome
INPUT_DIR=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2
IN_QUERY=($(find $INPUT_DIR -maxdepth 1 -name "mm10.fa"))
FILE_QUERY=${IN_QUERY[0]}
BASE_QUERY=${FILE_QUERY#${INPUT_DIR}/}
BASE_QUERY=${BASE_QUERY%.fa}

# query .gtf
INPUT_DIR=/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.99.GRCm38.p6.20200415.UCSCseqnames.gtf
IN_GTF=($(find $INPUT_DIR -maxdepth 1 -name "*gtf"))
FILE_GTF=${IN_GTF[0]}
BASE_GTF=${FILE_GTF#${INPUT_DIR}/}
BASE_GTF=${BASE_GTF%.gtf}

# script
SCRIPT=/common/WORK/fhorvat/programi/python/packages/bin/liftoff

# ----------------Commands------------------- #
# extract introns and coverage from bam file
$SCRIPT \
-g ${FILE_GTF} \
-p ${THREADS} \
-o ./${BASE_GTF}.${BASE_QUERY}.liftoff.gff \
${FILE_GENOME} \
${FILE_QUERY}
