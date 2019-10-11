#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/Dicer_Mili_KO/merged_replicates
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/datasets/2019_Sep/Data/Mapped/STAR_mm10/4_merged_replicates

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
