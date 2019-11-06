#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=~/public_html/Svoboda/bw_tracks/in_house/hamster_KO/hamster_testis_Mov10l.small_RNAseq/perfect_alignments.all_multimappers
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Mapped/bbmap_mesAur1

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find ${LOBSANG_PATH}/6_tracks -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL | sed '1d' > log.tracks_URL.txt
sed -i '1d' log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
