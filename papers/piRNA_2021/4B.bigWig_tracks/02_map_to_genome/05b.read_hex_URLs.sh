#/bin/bash
# ----------------Loading variables------------------- #
# get table URL
TABLE_PATH=http://hex.bioinfo.hr/~fhorvat/${HEX_PATH#~/public_html}/log.tracks_URL.txt

# ----------------Commands------------------- #
# download table from hex
wget ${TABLE_PATH}

# change permissions
find . -maxdepth 1 \( -name "*.bw" \) -exec chmod 744 {} \;
find . -maxdepth 1 \( -name "*.bam*" -not -name "*genome*" -not -name "*rDNA_45S*" -not -name "*sortedByName*" \) -exec chmod 744 {} \;
