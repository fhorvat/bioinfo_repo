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
INPUT_DIR=/common/DB/genome_reference/cow/bosTau9.ARS-UCD1.2.GCA_002263795.2
IN_GENOME=($(find $INPUT_DIR -maxdepth 1 -name "bosTau9.fa"))
FILE_GENOME=${IN_GENOME[0]}
BASE_GENOME=${FILE_GENOME#${INPUT_DIR}/}
BASE_GENOME=${BASE_GENOME%.fa}

# query genome
INPUT_DIR=/common/DB/genome_reference/cow/Btau_5.0.1.GCA_000003205.6
IN_QUERY=($(find $INPUT_DIR -maxdepth 1 -name "Btau_5.0.1.fa"))
FILE_QUERY=${IN_QUERY[0]}
BASE_QUERY=${FILE_QUERY#${INPUT_DIR}/}
BASE_QUERY=${BASE_QUERY%.fa}

# query .gtf
INPUT_DIR=.
IN_GTF=($(find $INPUT_DIR -maxdepth 1 -name "miRBase.22.Btau_5.0.1.20200202.genBankseqnames.gff3"))
FILE_GTF=${IN_GTF[0]}
BASE_GTF=${FILE_GTF#${INPUT_DIR}/}
BASE_GTF=${BASE_GTF%.gff3}

# script
SCRIPT=/common/WORK/fhorvat/programi/python/packages/bin/liftoff

# ----------------Commands------------------- #
# rename attributes in .gtf => "Derives_from" becomes "Parent"
sed -e 's/Derives_from/Parent/g' ${IN_GTF} > ${BASE_GTF}.fix.gff3

# extract introns and coverage from bam file
$SCRIPT \
-g ./${BASE_GTF}.fix.gff3 \
-p ${THREADS} \
-o ./${BASE_GTF}.${BASE_GENOME}.liftoff.gff \
-f ./feature_list.txt \
${FILE_GENOME} \
${FILE_QUERY}
