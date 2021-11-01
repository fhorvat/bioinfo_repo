#!/bin/bash
# ----------------Loading variables------------------- #
IN_DIR=.
OUT_DIR=.
FILE=($IN_DIR/knownGene*.txt.gz)
BASE=${FILE#$IN_DIR/}
BASE=${BASE%.txt.gz}
SCRIPT=/common/WORK/fhorvat/programi/UCSC/genePredToGtf

# ----------------Commands------------------- #
# converts genePred table from USCS's ftp server to gtf
zcat ${FILE} | cut -f1-10 | $SCRIPT file stdin stdout | sort -k1,1 -k4,4 | gzip > ${BASE}.gtf.gz
