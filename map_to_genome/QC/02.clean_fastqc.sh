#!/bin/bash

# ----------------Loading variables------------------- #
INPUT_DIR=.
EXP_NAME=`basename ${PWD%/Data/Raw/QC}`

# ----------------Commands------------------- #
# multiqc
multiqc .

# add to archive, remove
zip ${EXP_NAME}.zip *_fastqc.html
rm $INPUT_DIR/s_*.zip
rm $INPUT_DIR/*_fastqc.html
rm $INPUT_DIR/pbs.0* 
