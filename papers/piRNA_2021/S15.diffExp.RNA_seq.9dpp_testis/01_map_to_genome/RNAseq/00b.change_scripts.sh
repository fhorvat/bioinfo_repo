#/bin/bash
# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh
 
# get number of samples for PBS array index
N_SAMPLES=($(printf "%s\n" "${IN_SEQ[@]%_[1-2].txt.gz}" | sort -u | wc -l | awk '{print $1-1}'))

# get timestamp
TIMESTAMP=$(date)

# print all inputs/outputs to log
echo -e \
Mapping to genome pipeline started at ${TIMESTAMP} with following inputs and outputs:'\n'\
'\n'\
EXPERIMENT=$EXPERIMENT'\n'\
EXPERIMENT_NAME=$EXPERIMENT_NAME'\n'\
SINGLE_END=$SINGLE_END'\n'\
STAR_INDEX=$STAR_INDEX'\n'\
N_SAMPLES=`echo "$N_SAMPLES" | awk '{print $0+1}'`'\n'\
INPUT_DIR=$INPUT_DIR'\n'\
IN_SEQ=`printf "%s " "${IN_SEQ[@]}"`'\n'\
FEATURES_EXONS=$FEATURES_EXONS'\n'\
FEATURES_RMSK=$FEATURES_RMSK'\n'\
HEX_PATH=$HEX_PATH'\n'\
LOBSANG_PATH=$LOBSANG_PATH'\n' > log.mapping_pipeline.txt

# ----------------Commands------------------- #
# change scripts in loop
for script in "${SCRIPTS[@]}"
do
        echo $script
        perl -pi -e "s|%N_SAMPLES|$N_SAMPLES|" $script
        perl -pi -e "s|%HEX_PATH|$HEX_PATH|" $script
        perl -pi -e "s|%LOBSANG_PATH|$LOBSANG_PATH|" $script
done
