#/bin/bash
# ----------------Loading variables------------------- #
### source
source ../0_scripts/000.load_variables.sh

#### override source
# input bam files
INPUT_DIR=..
IN_SEQ=($(find ${INPUT_DIR} -maxdepth 1 \( -name "*.bam" \)))
REPLICATES=(${IN_SEQ[@]%_r[1-9]*.[S,P]E.bam})
REPLICATES=(${REPLICATES[@]#${INPUT_DIR}/})
REPLICATES=($(printf "%s\n" ${REPLICATES[@]} | sort -u))
N_REPLICATES=$((${#REPLICATES[@]} - 1))

# hex out path
HEX_OUT=accessory_data_sets/${EXPERIMENT}/merged_replicates
HEX_PATH=~/public_html/Svoboda/bw_tracks/${HEX_OUT}
LOBSANG_PATH=${PWD/common/common-lobsang}

### load
# scripts to change
SCRIPTS=($(find . -maxdepth 1 \( -name "0*.sh" -not -name "00b.change_scripts.sh" \)))
