THREADS=4
MEMORY=1g

ADAPTERS="/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Raw/mouse_B6.mm10/Cleaned/trimmed/3_logs/s_mouse_B6_r1.PE.adapters.fasta"

BBDUK_PAR="overwrite=t \
ktrim=r \
k=23 \
rcomp=t \
mink=12 \
hdist=1 \
minoverlap=8 \
minlength=25 \
tbo \
threads=$THREADS"

trim_PE = {
   
  // trim paired-end reads using bbduk
  doc "Trim poor quility bases and/or adapter sequences from reads"

  filter("trim","trim") {
    from("txt.gz","txt.gz") {
      exec """
        bbduk.sh
        in1=$input1.gz in2=$input2.gz
        out1=$output1.gz out2=$output2.gz
        ref=$ADAPTERS
        stats=${prefix}.stats 
        $BBDUK_PAR 
      """
     }
  }
}

run {
  "%_*.txt.gz" * [ trim_PE ]
}
