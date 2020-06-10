#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=40g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.01.correct_kmers.SE
#PBS -j oe
#PBS -J 0-1
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
THREADS=6
MEMORY=40g

INPUT_DIR=../../Raw/Cleaned/merged
IN_SEQ=(`ls ${INPUT_DIR}/*.txt.gz`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE/txt.gz/fastq}

# ----------------Commands------------------- #
# unzip the files
unpigz -p $THREADS -c $FILE > $BASE 

# correct erroneous k-mers
perl ~/bin/run_rcorrector.pl -s ${BASE} -t $THREADS
BASE2=${BASE/fastq/cor.fq}

# remove reads with "unfixable error" in header, strip "cor" from headers
sed -i -e '/unfixable_error/,+3d; s/ cor//g' ${BASE2}

# zip
pigz -p $THREADS ${BASE2}

# rename
mv ${BASE2}.gz ${BASE2/cor.fq/txt}.gz

# remove original unzipped file
[ -f "${BASE/fastq/txt.gz}" ] && rm ${BASE}
