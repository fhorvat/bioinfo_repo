#!/bin/bash

# ----------------Loading variables------------------- #
INPUT_DIR=.
EXP_NAME=`basename ${PWD%/Data/Raw/QC}`

# ----------------Commands------------------- #
# multiqc
multiqc .

# add to archive, remove
tar zcvf ${EXP_NAME}.fastqc.tar.gz *_fastqc.html
rm $INPUT_DIR/*.zip
rm $INPUT_DIR/*_fastqc.html

