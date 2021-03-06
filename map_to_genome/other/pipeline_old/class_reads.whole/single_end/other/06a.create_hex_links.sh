#/bin/bash

# ----------------Loading variables------------------- #
HEX_PATH=%HEX_PATH
LOBSANG_PATH=%LOBSANG_PATH

# ----------------Commands------------------- #
# create links
mkdir -p $HEX_PATH
cd $HEX_PATH
find $LOBSANG_PATH -name "*genome*.bw" -exec ln -s {} $HEX_PATH \;
find $LOBSANG_PATH -name "*genome*.bam*" -exec ln -s {} $HEX_PATH \;

# getURL
~/bin/getURL | sed '1d' > log.tracks_URL.txt
sed -i '1d' log.tracks_URL.txt
chmod 744 log.tracks_URL.txt

