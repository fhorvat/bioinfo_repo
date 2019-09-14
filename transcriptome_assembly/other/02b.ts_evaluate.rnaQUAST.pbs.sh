#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -l select=ncpus=6:mem=30g
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.02.rnaquast_ts_evaluate
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# script options
THREADS=6
RNA_QUAST=/common/WORK/fhorvat/programi/rnaQUAST/rnaQUAST-1.5.1/rnaQUAST.py

# fasta files
INPUT_DIR=.
IN_SEQ=(`ls ${INPUT_DIR}/*.fasta`)
FILES=$(printf " %s" ${IN_SEQ[@]})

# raw reads
BASE=${PWD/trinity/}
BASE=${BASE/spades/}
BASE=`basename $BASE`
INPUT_DIR=../../../Raw/$BASE/Cleaned
LEFT_READS=(${INPUT_DIR}/*all.PE_1.txt.gz)
RIGHT_READS=(${INPUT_DIR}/*all.PE_2.txt.gz)

# genome files
INPUT_DIR=../../../genomes/$BASE
GENOME_FASTA=(`find -L $INPUT_DIR -name "*.fa.gz" -or -name "*.fa"`)
GENOME_GTF=(`find -L $INPUT_DIR \( -name "ensembl.93*.gtf.gz" -or -name "ensembl.93*gtf" \) -and -name "*UCSCseqnames.gtf*"`)

if [ -z $GENOME_GTF ]
then
    GENOME_GTF=(`find -L $INPUT_DIR -name "ensembl.93*.gtf.gz" -or -name "ensembl.93*gtf"`)
fi

# ----------------Commands------------------- #
# ungzip fasta and gtf
if [[ $GENOME_FASTA =~ \.gz$ ]]
then
   unpigz -p $THREADS $GENOME_FASTA
   GENOME_FASTA=${GENOME_FASTA%.gz}
fi

if [[ $GENOME_GTF =~ \.gz$ ]]
then
   unpigz -p $THREADS $GENOME_GTF
   GENOME_GTF=${GENOME_GTF%.gz}
fi

# evaluate transcriptome assembly fasta files with rnaQUAST
python $RNA_QUAST \
--reference $GENOME_FASTA \
--gtf $GENOME_GTF \
--disable_infer_genes \
--disable_infer_transcripts \
--transcripts $FILES \
--left_reads $LEFT_READS \
--right_reads $RIGHT_READS \
--threads $THREADS \
--output_dir ./rnaQUAST_results

# gzip fasta and gtf
if [[ $GENOME_FASTA =~ \.fa$ ]]
then
   pigz -p $THREADS $GENOME_FASTA
fi

if [[ $GENOME_GTF =~ \.gtf$ ]]
then
   pigz -p $THREADS $GENOME_GTF
fi
