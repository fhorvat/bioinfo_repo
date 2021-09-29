#/bin/bash
# ----------------Loading variables------------------- #
### main dir
MAIN_DIR=../..

### source
source ${MAIN_DIR}/0_scripts/000.load_variables.sh

#### override source
# input bam files
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 \( -name "*.bam" \)))
REPLICATES=(${IN_SEQ[@]%_r[0-9]*.[S,P]E*.bam})
REPLICATES=(${REPLICATES[@]#${INPUT_DIR}/})
REPLICATES=($(printf "%s\n" ${REPLICATES[@]} | sort -u))
N_REPLICATES=$((${#REPLICATES[@]}))

# log files
LOG=$(find ${MAIN_DIR}/3_logs -maxdepth 1 -name "log.read_stats.txt")

# read lengths to loop through
READ_LENGTHS=(19to32nt 21to23nt 24to31nt all)
N_READ_LENGTHS=$((${#READ_LENGTHS[@]}))

# calculate number of jobs for PBS_ARRAY_INDEX
N_JOBS=$((${N_REPLICATES} * ${N_READ_LENGTHS} - 1))
 
# get layout (single or paired end)
if [ ${SINGLE_END} == TRUE ]
then
        LAYOUT="SE"
else
        LAYOUT="PE"
fi

# hex out path
HEX_OUT=${HEX_OUT}/merged_replicates
HEX_PATH=~/public_html/Svoboda/bw_tracks/${HEX_OUT}
LOBSANG_PATH=${PWD/common/common-lobsang}

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))
