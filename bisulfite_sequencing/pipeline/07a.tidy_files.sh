#/bin/bash
# ----------------Commands------------------- #
### create dirs
mkdir \
0_scripts \
1_mapped \
2_unmapped \
3_logs \
4_methylation_calls \
5_reports

### move files
# original scripts
mv \
0*.sh \
0*.R \
0*.txt \
./0_scripts/

# bams
#find . -maxdepth 1 \( -name "*genome*.bam" -or -name "*45S*.bam" -or -name "*sortedByName.bam" \) -exec mv {} 1_mapped/ \;

# unmapped files
mv \
*.unmapped.txt.gz \
./2_unmapped/

# other logs and scripts
mv \
log.tracks_URL.txt \
log.mapping_pipeline.txt \
pbs*.* \
./3_logs/

# methylation calls
mv \
CHG_context_* \
CHH_context_* \
CpG_context_* \
./4_methylation_calls/

# reports
mv \
*_bismark_bt2_*_report.txt \
*.bismark.cov \
*.M-bias.txt \
*_splitting_report.txt \
bismark_summary_report.txt \
./5_reports/
