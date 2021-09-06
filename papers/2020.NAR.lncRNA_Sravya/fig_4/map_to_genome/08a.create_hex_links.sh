#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/maternal_transcriptomes/rat.rn6
LOBSANG_PATH=/common-lobsang/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/rat.rn6

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL | sed '1d' > log.tracks_URL.txt
sed -i '1d' log.tracks_URL.txt
chmod 744 log.tracks_URL.txt
