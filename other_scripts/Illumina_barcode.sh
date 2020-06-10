zgrep "^@" $sample | cut -d":" -f10 | sort -r | uniq -c | sort -nrk1,1 > $sample.barcode
