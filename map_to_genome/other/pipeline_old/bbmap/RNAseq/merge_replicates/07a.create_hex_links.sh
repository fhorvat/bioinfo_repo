#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=%HEX_PATH/merged_replicates
LOBSANG_PATH=%LOBSANG_PATH/4_merged_replicates

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
