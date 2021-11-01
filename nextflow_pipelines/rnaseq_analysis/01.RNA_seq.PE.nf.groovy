#!/usr/bin/env nextflow

//// variables and scripts
INPUT_DIR = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/test"
READS = "$INPUT_DIR/s_*{1,2}.txt.gz"

GENOME_DIR = "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
SJDB_OVERHANG = "sjdbOverhang_249"
STAR_INDEX = "$GENOME_DIR/STAR_index/$SJDB_OVERHANG"
CHR_LENGTH = "$STAR_INDEX/chrNameLength.txt"
GTF_PATH = "$GENOME_DIR/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz"
RMSK_PATH = "$GENOME_DIR/rmsk.mm10.20180919.clean.fa.out.gz"
GENES_INFO_PATH = "$GENOME_DIR/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv"
SAMPLE_TABLE_PATH = "$INPUT_DIR/../../Documentation/CNOT6L.sample_table.csv"
GROUPING_VARIABLES = "stage genotype"
RESULTS_GROUPS = "GV_KO,GV_WT MII_KO,MII_WT 1C_KO,1C_WT"

HEX_PATH = "~/public_html/Svoboda/bw_tracks/accessory_data_sets/nextflow_test"
LOBSANG_PATH = "${PWD}".replaceFirst(/common/, "common-lobsang")
TABLE_PATH = "http://hex.bioinfo.hr/~fhorvat" + "${HEX_PATH}".replaceFirst(/~\/public_html/, "") + "/log.tracks_URL.txt"

SCRIPTS_DIR = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/test/bin"
CLASS_READS_SCRIPT = "$SCRIPTS_DIR/class_reads.R"
SUM_READ_STATS_SCRIPT = "$SCRIPTS_DIR/sum_read_stats.R"
SCALE_BIGWIG_SCRIPT= "$SCRIPTS_DIR/scale_bigWig.R"
STATS_AND_TRACKS_SCRIPT = "$SCRIPTS_DIR/write_stats_and_tracks.R"
SUMMARIZE_OVERLAPS_SCRIPT = "$SCRIPTS_DIR/summarizeOverlaps.R"
FPKM_TABLES_SCRIPT = "$SCRIPTS_DIR/FPKM_tables.R"
DIFFEXP_ANALYSIS_SCRIPT = "$SCRIPTS_DIR/diffExp_analysis.R"

//// channels
// fastq to fastqc channel
Channel
  .fromFilePairs( "${READS}" )
  .ifEmpty { error "Cannot find any reads matching: ${READS}" }
  .set { read_pairs_for_fastqc }


// fastq to mapping channel
Channel
  .fromFilePairs( "${READS}", flat: true )
  .ifEmpty { error "Cannot find any reads matching: ${READS}" }
  .set { read_pairs_for_map }


//// processes
// map reads to the genome
process map_to_genome{

  publishDir "results.mapped_reads", mode: "copy", pattern: "{*bam,*bam.bai}"

  input:
    set val(sample_id), file(reads_1), file(reads_2) from read_pairs_for_map

  output:
    set val(sample_id), file("${sample_id}.bam") into mapped_bams_for_tracks, mapped_bams_for_sort, mapped_bams_for_bamtools_stats
    set val(sample_id), file("${sample_id}.bam.bai") into bam_indices
    file "*.Log.final.out" into star_logs_for_stats, star_logs_for_multiqc
    file "*.ReadsPerGene.out.tab" into star_counts_for_multiqc
    file "*.bam" into mapped_bams_for_summarizedOverlaps

  // PBS parameters
  cpus 2
  memory "40 GB"
  executor "pbspro"

  // run
  """
  # map to genome
  STAR \\
  --genomeDir $STAR_INDEX \\
  --readFilesIn $reads_1 $reads_2 \\
  --runThreadN 2 \\
  --genomeLoad NoSharedMemory \\
  --limitBAMsortRAM  20000000000 \\
  --readFilesCommand unpigz -c \\
  --outFileNamePrefix ${sample_id}. \\
  --outSAMtype BAM SortedByCoordinate \\
  --outReadsUnmapped Fastx \\
  --outFilterMultimapNmax 99999 \\
  --outFilterMultimapScoreRange 0 \\
  --outFilterMismatchNoverLmax 0.05 \\
  --sjdbScore 2 \\
  --quantMode GeneCounts

  # rename
  mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam

  # index files
  samtools index ${sample_id}.bam
  """

}


// creates raw bigWig tracks
process bigwig_tracks{

  publishDir "results.tracks_raw", mode: "copy"

  input:
    set val(sample_id), file(mapped_bam) from mapped_bams_for_tracks

  output:
    set val(sample_id), file("${sample_id}.bw") into raw_tracks_for_scale

  // PBS parameters
  cpus 2
  memory "10 GB"
  executor "pbspro"

  // run
  """
  # produce bigWig tracks
  genomeCoverageBed -ibam ${sample_id}.bam -bg -split -g $CHR_LENGTH > ${sample_id}.bedGraph
  wigToBigWig ${sample_id}.bedGraph $CHR_LENGTH ${sample_id}.bw
  [ -f "${sample_id}.bw" ] && rm ${sample_id}.bedGraph
  """

}


// sort joined bam by name
process sort_bam{

  input:
    set val(sample_id), file(mapped_bam) from mapped_bams_for_sort

  output:
    set val(sample_id), file("${sample_id}.sortedByName.bam") into sorted_bams

  // PBS parameters
  cpus 2
  memory "10 GB"
  executor "pbspro"

  // run
  """
  # sort bam by name
  sambamba sort -m 10GB -o ${sample_id}.sortedByName.bam -n -t 2 ${sample_id}.bam
  """

}


// class reads in sorted bams
process class_reads{

  publishDir "results.stats_and_tracks", mode: "copy"

  input:
    set val(sample_id), file(sorted_bam) from sorted_bams

  output:
    file "*.read_stats.txt" into read_stats_for_summary
    set val(sample_id), file("${sample_id}.read_stats.txt") into read_stats_for_scale

  // PBS parameters
  cpus 1
  memory "20 GB"
  executor "pbspro"

  // run
  """
  # class reads
  Rscript $CLASS_READS_SCRIPT \\
  --bam_file ${sample_id}.sortedByName.bam \\
  --gtf_path ${GTF_PATH} \\
  --rmsk_path ${RMSK_PATH}
  """

}


// summarizes stats of reads in all .bams
process sum_read_stats{

  publishDir "results.stats_and_tracks", mode: "copy"

  input:
    file(read_stats_list) from read_stats_for_summary.collect()
    file(star_logs_list) from star_logs_for_stats.collect()

  output:
    file "log.read_stats.txt" into read_stats_sum

  // PBS parameters
  //cpus 1
  //memory "1 GB"
  //executor "pbspro"

  // run
  """
  # summarize read class stats
  Rscript $SUM_READ_STATS_SCRIPT \\
  --read_stats_list ${(read_stats_list)} \\
  --mapped_logs_list ${(star_logs_list)}
  """

}


// combine tracks and stats to one channel
Channel
  raw_tracks_for_scale.join(read_stats_for_scale)
  .set { stats_and_tracks }

// scale bigWig to RPMs
process scale_bigWig{

  publishDir "results.tracks_scaled", mode: "copy"

  input:
    set val(sample_id), file("${sample_id}.bw"), file("${sample_id}.read_stats.txt") from stats_and_tracks

  output:
    file "*.scaled.bw" into scaled_tracks_for_links

  // PBS parameters
  cpus 1
  memory "10 GB"
  executor "pbspro"

  // run
  """
  Rscript $SCALE_BIGWIG_SCRIPT \\
  --raw_tracks ${sample_id}.bw \\
  --read_stats ${sample_id}.read_stats.txt
  """

}


// make links on hex.bioinfo.hr
process make_links{

  publishDir "results.stats_and_tracks", mode: "copy"

  input:
    file scaled_track from scaled_tracks_for_links.collect()

  output:
    file "log.tracks_URL.txt" into tracks_URLs

  // run
  """
  # create links on hex
  ssh -t fhorvat@hex.bioinfo.hr 'mkdir -p $HEX_PATH; \\
  cd $HEX_PATH; \\
  ln -s $LOBSANG_PATH/results.mapped_reads/* .; \\
  ln -s $LOBSANG_PATH/results.tracks_raw/* .; \\
  ln -s $LOBSANG_PATH/results.tracks_scaled/* .; \\
  ~/bin/getURL > log.tracks_URL.txt; \\
  chmod 744 log.tracks_URL.txt'

  # download links from hex
  wget $TABLE_PATH

  # change permissions
  find ${PWD}/results.tracks_raw -name "*.bw" -exec chmod 744 {} +
  find ${PWD}/results.tracks_scaled -name "*.bw" -exec chmod 744 {} +
  find ${PWD}/results.mapped_reads -name "*.bam*" -exec chmod 744 {} +
  """

}


// create stats and tracks table
process write_stats_and_tracks{

  publishDir "results.stats_and_tracks", mode: "copy"

  input:
    file read_stats_sum from read_stats_sum
    file tracks_URLs from tracks_URLs

  output:
    file "*.stats_and_tracks.csv" into final_stats_and_tracks

  // run
  """
  Rscript $STATS_AND_TRACKS_SCRIPT \\
  --read_stats_sum ${read_stats_sum} \\
  --tracks_URLs ${tracks_URLs}
  """

}


// does quality check on fastq files for multiQC
process fastqc{

  input:
    set val(sample_id), file(reads) from read_pairs_for_fastqc

  output:
    file "*_fastqc.{zip,html}" into fastqc_results_for_multiqc

  // PBS parameters
  cpus 1
  memory "1 GB"
  executor "pbspro"

  script:
  """
  fastqc -q $reads
  """

}


// creates bamtools stats report for multiQC
process bamtools_stats{

  input:
    set val(sample_id), file(mapped_bam) from mapped_bams_for_bamtools_stats

  output:
    file "*stats" into bamtools_stats_for_multiqc

  // PBS parameters
  cpus 1
  memory "5 GB"
  executor "pbspro"

  // run
  """
  bamtools stats -in ${sample_id}.bam | tail -n +2 > ${sample_id}.stats
  """

}


// creates QC report
process multiqc{

  publishDir "results.multiqc", mode: "copy"

  input:
    file fastqc_results from fastqc_results_for_multiqc.collect().ifEmpty([])
    file star_logs from star_logs_for_multiqc.collect().ifEmpty([])
    file star_counts from star_counts_for_multiqc.collect().ifEmpty([])
    file bamtools_stats from bamtools_stats_for_multiqc.collect().ifEmpty([])

  output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

  script:
  """
  multiqc .
  """

}


// summarizeOverlaps with Ensembl .gtf
process summarizeOverlaps{

  input:
    file mapped_bams_list from mapped_bams_for_summarizedOverlaps.collect()

  output:
    file "*.se.RDS" into summarizeOverlaps_for_FPKM, summarizeOverlaps_for_diffExp

  // PBS parameters
  cpus 1
  memory "10 GB"
  executor "pbspro"

  // run
  """
  # summarize read class stats
  Rscript $SUMMARIZE_OVERLAPS_SCRIPT \\
  --mapped_bams_list ${(mapped_bams_list)} \\
  --gtf_path ${GTF_PATH}
  """

}


// calculate FPKM and average FPKM
process FPKM_tables{

  publishDir "results.FPKM_tables", mode: "copy"

  input:
    file se_rds from summarizeOverlaps_for_FPKM

  output:
    file "*.FPKM.csv"
    file "*.FPKM_long.csv"
    file "*.FPKM_mean.csv" into FPKM_table_for_diffExp
    file "*.FPKM_stats.csv"

  // PBS parameters
  cpus 1
  memory "10 GB"
  executor "pbspro"

  // run
  """
  # summarize read class stats
  Rscript $FPKM_TABLES_SCRIPT \\
  --se_path ${se_rds} \\
  --sample_table_path ${SAMPLE_TABLE_PATH} \\
  --gtf_path ${GTF_PATH} \\
  --genes_info_path ${GENES_INFO_PATH} \\
  --grouping_variables ${GROUPING_VARIABLES}
  """

}


// does differential expression analysis along with some exploratory plots (PCA, heatmaps, etc.)does differential expression analysis along with some exploratory plots (PCA, heatmaps, etc.)
//
process diffExp_analysis{

  publishDir "results.diffExp", mode: "copy"

  input:
    file se_rds from summarizeOverlaps_for_diffExp
    file fpkm_tb from FPKM_table_for_diffExp

  output:
    file "*.PCA_plot.png"
    file "*.PCA_plot.html"
    file "*.distance_heatmap.png"
    file "*.DESeq2.all_results.xlsx"
    file "*.DESeq2.significant_results.xlsx"
    file "*.MA_plot.png"
    file "*.MA_plot.html"
    file "*.vulcano_plot.png"
    file "*.vulcano_plot.html"

  // PBS parameters
  cpus 1
  memory "10 GB"
  executor "pbspro"

  // run
  """
  Rscript $DIFFEXP_ANALYSIS_SCRIPT \\
  --se_path ${se_rds} \\
  --fpkm_path ${fpkm_tb} \\
  --sample_table_path ${SAMPLE_TABLE_PATH} \\
  --genes_info_path ${GENES_INFO_PATH} \\
  --grouping_variables ${GROUPING_VARIABLES} \\
  --results_groups ${RESULTS_GROUPS}
  """

}

