#/bin/bash
# ----------------Loading variables------------------- #
### main dir
MAIN_DIR=..

### source
source ${MAIN_DIR}/0_scripts/000.load_variables.sh

#### override source
# input bam files
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 \( -name "*.bam" \)))
REPLICATES=(${IN_SEQ[@]%_r[1-9]*.[S,P]E.bam})
REPLICATES=(${REPLICATES[@]#${INPUT_DIR}/})
REPLICATES=($(printf "%s\n" ${REPLICATES[@]} | sort -u))
N_JOBS=$((${#REPLICATES[@]} - 1))

# get layout (single or paired end)
if [ ${SINGLE_END} == TRUE ]
then
        LAYOUT="SE"
else
        LAYOUT="PE"
fi

# log file
LOG=${MAIN_DIR}/3_logs/log.read_stats.txt

# hex out path
HEX_OUT=${HEX_OUT}/merged_replicates
HEX_PATH=~/public_html/Svoboda/bw_tracks/${HEX_OUT}
LOBSANG_PATH=${PWD/common/common-lobsang}

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))
