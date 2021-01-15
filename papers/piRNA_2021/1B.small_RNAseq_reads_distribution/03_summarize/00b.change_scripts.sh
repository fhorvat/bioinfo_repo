#/bin/bash
# ----------------Loading variables------------------- #
# source
source ./000.load_variables.sh
 
# get number of samples for PBS array index
N_SAMPLES=($(printf "%s\n" "${IN_SEQ[@]%bam}" | sort -u | wc -l | awk '{print $1-1}'))

# get timestamp
TIMESTAMP=$(date)

# ----------------Commands------------------- #
# change scripts in loop
for script in "${SCRIPTS[@]}"
do
        echo $script
        perl -pi -e "s|%N_SAMPLES|$N_SAMPLES|" $script
done
