// ----------------Load variables-------------------
GENOME_DIR="/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"
SJDB_OVERHANG="sjdbOverhang_249"
STAR_INDEX="$GENOME_DIR/STAR_index/$SJDB_OVERHANG"
CHR_LENGTH="${STAR_INDEX}/chrNameLength.txt"

// ----------------Define runs-------------------
align_STAR = {
  produce("*.Aligned.sortedByCoord.out.bam"){
    exec """
      STAR
      --genomeDir $STAR_INDEX
      --readFilesIn $input1.gz $input2.gz
      --runThreadN $threads
      --genomeLoad LoadAndRemove
      --limitBAMsortRAM  20000000000
      --readFilesCommand unpigz -c
      --outFileNamePrefix $output.prefix.
      --outSAMtype BAM SortedByCoordinate
      --outReadsUnmapped Fastx
      --outFilterMultimapNmax 99999
      --outFilterMultimapScoreRange 0
      --outFilterMismatchNoverLmax 0.2
      --sjdbScore 2
    """
  }
}

// ----------------Run pipeline-------------------
// s_GV_WT_r1.PE_1.txt.gz

run { 
  "s_%.PE_*.txt.gz" * [ align_STAR ] 
}
