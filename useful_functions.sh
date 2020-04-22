#!/bin/bash
# convert geneBank format to fasta
gb2fasta () {
   awk '/^LOCUS   / {printf(">%s\n",$2);next;} /^ORIGIN/ {inseq=1;next;} /^\/\// {inseq=0;} {if(inseq==0) next; gsub(/[0-9 ]/,"",$0); printf("%s\n",toupper($0));}' $1
}

# convert blast output 6 to bed 
blast2bed () {
   awk '{if ($9 < $10) print $2"\t"$9-1"\t"$10"\t"$1"\t"0"\t+"; else print $2"\t"$10-1"\t"$9"\t"$1"\t"0"\t-"}' $1
}

# get bed from fasta
fasta2bed () {
   samtools faidx $1
   awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${1}.fai > ${1/fasta/bed}
}

# remove line breaks from fasta
removeFastaBreaks () {
   awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' $1
}

# zcat + head
zhead () {
   zcat $1 | head
}

# get fasta from fastq
#fastq2fasta () {
#   zcat $1 | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n"
#}

# get fasta from fastq using sed
fastq2fasta () {
   zcat $1 | sed -n '1~4s/^@/>/p;2~4p'
}

# count reads in fastq
countFastq () {
   zcat $1 | awk '{s++}END{print s/4}'
}

# get length of variable (bash array)
vlen () {
    declare ARRAY_NAME="$1"
    declare INDIRECT_REFERENCE="\${#${ARRAY_NAME}[@]}"
    case "$-" in
    *'u'*)
        set +u
        eval "echo \"${INDIRECT_REFERENCE}\""
        set -u
        ;;
    *)
        eval "echo \"${INDIRECT_REFERENCE}\""
        ;;
    esac
}

# print whole array
vprnt () {
    declare ARRAY_NAME="$1"
    declare INDIRECT_REFERENCE="\${${ARRAY_NAME}[@]}"
    case "$-" in
    *'u'*)
        set +u
        eval "echo \"${INDIRECT_REFERENCE}\""
        set -u
        ;;
    *)
        eval "echo \"${INDIRECT_REFERENCE}\""
        ;;
    esac
}

# prints csv in tidy way
tcsv () {
   column -s, -t < $1 | head
}
