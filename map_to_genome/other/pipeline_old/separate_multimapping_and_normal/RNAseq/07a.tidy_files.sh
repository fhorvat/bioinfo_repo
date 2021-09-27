#/bin/bash
# ----------------Commands------------------- #
### create dirs
mkdir \
0_scripts \
1_mapped \
2_unmapped \
3_logs

### move files
# original scripts
mv \
0*.sh \
0*.R \
0*.txt \
0_scripts/

# bams
find . -maxdepth 1 \( -name "*genome*.bam" -or -name "*45S*.bam" -or -name "*sortedByName.bam" \) -exec mv {} 1_mapped/ \;

# unmapped files
mv \
*Unmapped.out* \
2_unmapped/

# other logs and scripts
mv \
*.read_stats.txt \
log.tracks_URL.txt \
log.mapping_pipeline.txt \
*.out \
*.tab \
pbs*.* \
3_logs/
