#/bin/bash
# ----------------Loading variables------------------- #
HEX_PATH=/home/students/fhorvat/public_html/Svoboda/bw_tracks/in_house/lncRNA_KO/Lnc1_KO.2018_Dec
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
wget $TABLE_PATH
chmod 744 *genome*bam* *genome*bw

