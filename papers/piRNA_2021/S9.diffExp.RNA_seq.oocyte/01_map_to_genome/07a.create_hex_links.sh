#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/hamster_KO/hamster_oocyte_Mov10l.RNAseq
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Mapped/STAR_mesAur1

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL > log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
