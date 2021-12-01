find /common/WORK/fhorvat/Projekti/Svoboda -not -path '*vfranke*' \( -name "*.R" \) -exec grep "^library(.*)$" {} \; | sed 's/"//g' | sort | uniq -c | sort -k1,1 -n -r > scripts_counts.txt
