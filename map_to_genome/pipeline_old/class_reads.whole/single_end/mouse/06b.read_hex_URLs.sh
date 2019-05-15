#/bin/bash
# ----------------Loading variables------------------- #
HEX_PATH=%HEX_PATH
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
wget $TABLE_PATH
chmod 744 *genome*bam* *genome*bw

