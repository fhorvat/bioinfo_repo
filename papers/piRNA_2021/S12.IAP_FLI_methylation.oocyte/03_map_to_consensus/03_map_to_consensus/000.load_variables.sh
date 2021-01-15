#/bin/bash
# ----------------Loading variables------------------- #
#### experiment
EXPERIMENT=${PWD%/Data*}
EXPERIMENT=${EXPERIMENT##/*/}
EXPERIMENT_NAME=${EXPERIMENT%_*_*}

#### input
SINGLE_END=FALSE
INPUT_DIR=../../../filtered_reads/IAPLTR3

# input fastq files
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "s_*txt.gz" \)))

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))

# input genome files
BISMARK_INDEX=../../../bismark_index/FLI_consensus.IAPLTR3
CHR_LENGTH=${BISMARK_INDEX}/chrNameLength.txt
