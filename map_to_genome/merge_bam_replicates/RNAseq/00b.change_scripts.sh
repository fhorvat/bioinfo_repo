#/bin/bash
# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh
 
# get timestamp
TIMESTAMP=$(date)

# print all inputs/outputs to log
echo -e \
Merging replicates pipeline started at ${TIMESTAMP} with following inputs and outputs:'\n'\
'\n'\
EXPERIMENT=${EXPERIMENT}'\n'\
EXPERIMENT_NAME=${EXPERIMENT_NAME}'\n'\
STAR_INDEX=${STAR_INDEX}'\n'\
N_REPLICATES=`echo $((${N_REPLICATES} + 1)) `'\n'\
INPUT_DIR=${INPUT_DIR}'\n'\
REPLICATES=`printf "%s " "${REPLICATES[@]}"`'\n'\
HEX_PATH=$HEX_PATH'\n'\
LOBSANG_PATH=$LOBSANG_PATH'\n' > log.merging_replicates_pipeline.txt

# ----------------Commands------------------- #
# change scripts in loop
for script in "${SCRIPTS[@]}"
do
        echo $script
        perl -pi -e "s|%N_REPLICATES|$N_REPLICATES|" $script
        perl -pi -e "s|%HEX_PATH|$HEX_PATH|" $script
        perl -pi -e "s|%LOBSANG_PATH|$LOBSANG_PATH|" $script
done
