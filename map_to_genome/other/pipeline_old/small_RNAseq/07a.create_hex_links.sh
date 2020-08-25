#/bin/bash

# ----------------Loading variables------------------- #
EXPERIMENT=Yang_2016_SciAdv_GSE83581
HEX_PATH=~/public_html/Svoboda/bw_tracks/accessory_data_sets/${EXPERIMENT}
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/${EXPERIMENT}/Data/Mapped/STAR_mm10.21_to_23nt.perfect

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
