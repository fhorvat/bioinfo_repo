#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=1:mem=10g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.gtf_to_genePred
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -name "*liftoff.gff"))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.gff}

INPUT_DIR=../../
CHROM_SIZES=($(find ${INPUT_DIR} -name "*chrom.sizes"))

# ----------------Commands------------------- #
# gtf to BigGenePred 
gtfToGenePred -infoOut=${BASE}.info.txt ${FILE} ${BASE}.genePred

# convert the genePred extended file to a pre-bigGenePred text file
genePredToBigGenePred ${BASE}.genePred stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${BASE}.bigGenePred.txt

# convert text bigGenePred to a binary indexed format
bedToBigBed -type=bed12+8 -extraIndex=name -tab -as=/common/WORK/fhorvat/programi/UCSC/bigGenePred.as ${BASE}.bigGenePred.txt ${CHROM_SIZES} ${BASE}.bb

# set permission
chmod 744 ${BASE}.bb
