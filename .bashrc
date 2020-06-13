# Sample .bashrc for SuSE Linux
# Copyright (c) SuSE GmbH Nuernberg

# There are 3 different types of shells in bash: the login shell, normal shell
# and interactive shell. Login shells read ~/.profile and interactive shells
# read ~/.bashrc; in our setup, /etc/profile sources ~/.bashrc - thus all
# settings made here will also take effect in a login shell.
#
# NOTE: It is recommended to make language settings in ~/.profile rather than
# here, since multilingual X sessions would not work properly if LANG is over-
# ridden in every subshell.

# Some applications read the EDITOR variable to determine your favourite text
# editor. So uncomment the line below and enter the editor of your choice :-)
#export EDITOR=/usr/bin/vim
#export EDITOR=/usr/bin/mcedit

# For some news readers it makes sense to specify the NEWSSERVER variable here
#export NEWSSERVER=your.news.server

# If you want to use a Palm device with Linux, uncomment the two lines below.
# For some (older) Palm Pilots, you might need to set a lower baud rate
# e.g. 57600 or 38400; lowest is 9600 (very slow!)
#
#export PILOTPORT=/dev/pilot
#export PILOTRATE=115200

test -s ~/.alias && . ~/.alias || true

# ----------------ALIASES------------------- #
alias ls="ls --color -h --group-directories-first"
alias LS="ls"
alias CD="cd" 
alias CAT="cat"
alias QSTAT="qstat -at -u fhorvat"
alias dtrx="python2 /common/WORK/fhorvat/programi/dtrx/dtrx-7.1/bin/dtrx"
alias prnt="printf '%s\n'"
alias pip_install="pip install --install-option="--prefix=/common/WORK/fhorvat/programi/python/packages""
alias pip3_install="pip3 install --install-option="--prefix=/common/WORK/fhorvat/programi/python/packages""
alias matrix="cmatrix -sb -u 7"
alias qsub_interactive="qsub -I -N pbs.interactive -M fihorvat@gmail.com -m n -q MASTER -l select=ncpus=1:mem=100g"
alias DU="du -h -d 1 | sort -h -r"
alias t1="tail -n 1"
alias xc="xclip"
alias hex="ssh hex.bioinfo.hr"
alias tpbs="tail -n 1 pbs.*"
alias L="l"
alias sview="samtools view"
alias bsn="basename $PWD"

# path variables
alias go="cd /common/WORK/fhorvat/Projekti/Svoboda"
alias home="cd /common/WORK/fhorvat"
alias scripts="cd /common/WORK/fhorvat/Projekti/Svoboda/scripts"
alias sra="cd /common/DB/SRA/sra/Svoboda/2017_download"
alias ref="cd /common/DB/genome_reference"


# ----------------PATHS AND VARIABLES------------------- #
# my variables
export SCRIPTS=/common/WORK/fhorvat/Projekti/Svoboda/scripts

# Java
export JAVA_HOME=/common/software/JDK8/jdk1.8.0_20
export PATH=$JAVA_HOME/bin:$PATH

# man
MANPATH=$MANPATH:$HOME/share/man

# MOSAIK
export MOSAIK_TMP=/common/WORK/fhorvat/tmp/MOSAIK

# SHRiMP2
export SHRIMP_FOLDER=/common/WORK/fhorvat/programi/SHRiMP2/SHRiMP_2_2_3

# python
export PYTHONPATH=/common/WORK/fhorvat/programi/python/packages:$PYTHONPATH
export PYTHONUSERBASE=/common/WORK/fhorvat/programi/python/packages
export PATH=/common/WORK/fhorvat/programi/python/packages/bin:$PATH
export PATH=/common/WORK/fhorvat/programi/python/packages/lib/python2.7/site-packages/bin:$PATH
export PATH=/common/WORK/fhorvat/programi/python/packages/lib/python3.6/site-packages/bin:$PATH
export PYTHONPATH=${PYTHONPATH}:/common/WORK/fhorvat/programi/mafTools

# boost
export LD_LIBRARY_PATH=/common/WORK/fhorvat/programi/boost/boost_1_66_0/lib

# R common library path
export R_LIBS_USER=/common/WORK/fhorvat/programi/R/libraries

# UCSC tools
export PATH=/common/WORK/fhorvat/programi/UCSC:$PATH

# udunits
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/common/WORK/fhorvat/programi/udunits2/udunits-2.2.26/lib 

# Trinity
export TRINITY_HOME=/common/WORK/fhorvat/programi/Trinity/trinityrnaseq-Trinity-v2.8.4

# busco
export PATH=$PATH:/common/WORK/fhorvat/programi/busco/busco-master/scripts/
export PATH=$PATH:/common/WORK/fhorvat/programi/hmmer/hmmer-3.1b2/bin/
export PATH=$PATH:/common/WORK/kristian/bin/EMBOSS/EMBOSS-6.6.0/emboss/transeq

# miniconda3
export PATH=$PATH:/common/WORK/fhorvat/programi/Miniconda3/bin

# bpipe
export PATH=$PATH:/common/WORK/fhorvat/programi/bpipe/bpipe-0.9.9.7/bin

# groovy
export GROOVY_HOME=/common/WORK/fhorvat/programi/groovy/groovy-2.5.7
export PATH=$PATH:$GROOVY_HOME/bin

# xclip
export PATH=$PATH:/home/students/fhorvat/local/bin

# mugsy
export MUGSY_INSTALL=/common/WORK/fhorvat/programi/mugsy/mugsy_x86-64-v1r2.3
export PATH=$PATH:$MUGSY_INSTALL
export PERL5LIB=$MUGSY_INSTALL

# HMMER
export PATH=$PATH:/common/WORK/fhorvat/programi/hmmer/hmmer-3.2.1/bin

# bedtools
export PATH=$PATH:/common/WORK/fhorvat/programi/bedtools/bedtools-2.29.2/bin

# Perl
#PATH="/home/students/fhorvat/perl5/bin${PATH:+:${PATH}}"; export PATH;
#PERL5LIB="/home/students/fhorvat/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
#PERL_LOCAL_LIB_ROOT="/home/students/fhorvat/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
#PERL_MB_OPT="--install_base \"/home/students/fhorvat/perl5\""; export PERL_MB_OPT;
#PERL_MM_OPT="INSTALL_BASE=/home/students/fhorvat/perl5"; export PERL_MM_OPT;

PERL5LIB="/common/WORK/fhorvat/programi/Perl/libraries/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
PERL5LIB="/common/WORK/fhorvat/programi/Perl/libraries/lib/perl5/site_perl${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;

# NCBI
export LD_LIBRARY_PATH=/common/WORK/fhorvat/programi/ncbi/ngs/lib64:$LD_LIBRARY_PATH
export NGS_LIBDIR=/common/WORK/fhorvat/programi/ncbi/ngs/lib64
export LD_LIBRARY_PATH=/common/WORK/fhorvat/programi/ncbi/ncbi-vdb/lib64:$LD_LIBRARY_PATH

# SPAdes
export PATH=$PATH:/common/WORK/fhorvat/programi/SPAdes/SPAdes-3.14.1-Linux/bin

# TBB
export TBBROOT="/common/WORK/fhorvat/programi/TBB/oneTBB" #
tbb_bin="/common/WORK/fhorvat/programi/TBB/oneTBB/build/linux_intel64_gcc_cc4.8_libc2.19_kernel3.12.28_release" #
if [ -z "$CPATH" ]; then #
    export CPATH="${TBBROOT}/include" #
else #
    export CPATH="${TBBROOT}/include:$CPATH" #
fi #
if [ -z "$LIBRARY_PATH" ]; then #
    export LIBRARY_PATH="${tbb_bin}" #
else #
    export LIBRARY_PATH="${tbb_bin}:$LIBRARY_PATH" #
fi #
if [ -z "$LD_LIBRARY_PATH" ]; then #
    export LD_LIBRARY_PATH="${tbb_bin}" #
else #
    export LD_LIBRARY_PATH="${tbb_bin}:$LD_LIBRARY_PATH" #
fi #
 #

# AUGUSTUS and BRAKER
export BAMTOOLS_PATH=/common/WORK/fhorvat/programi/bamtools/bamtools/bin
export AUGUSTUS_CONFIG_PATH=/common/WORK/fhorvat/programi/Augustus/augustus-3.3.3/config
export AUGUSTUS_BIN_PATH=/common/WORK/fhorvat/programi/Augustus/augustus-3.3.3/bin
export AUGUSTUS_SCRIPTS_PATH=/common/WORK/fhorvat/programi/Augustus/augustus-3.3.3/scripts
export BLAST_PATH=/common/WORK/fhorvat/programi/BLAST+/ncbi-blast-2.9.0+/bin
export SAMTOOLS_PATH=/common/WORK/fhorvat/programi/samtools/samtools-1.3.1/bin
export GENEMARK_PATH=/common/WORK/jaruizs/programs/gm_et_linux_64/gmes_petap
export DIAMOND_PATH=/common/WORK/jaruizs/programs/diamond
export ALIGNMENT_TOOL_PATH=/common/WORK/jaruizs/programs/gth-1.7.1-Linux_x86_64-64bit/bin

# ----------------FUNCTIONS------------------- #
# load useful functions
source /common/WORK/fhorvat/Projekti/Svoboda/scripts/useful_functions.sh

# changes dir to n-th dir in current directory 
function cdd () { 
   cd $(ls -d */ | head -n $1 | tail -n 1)
}

# path function
function path(){
    old=$IFS
    IFS=:
    printf "%s\n" $PATH
    IFS=$old
}

# prompt
prompt_cmd() {
 [ -f ~/.x11 ] && . ~/.x11
 echo -ne "\033k${PWD##*/}\033\\"
}

# time stamp 
timestamp() {
  date +"%F_%T_%s"
}

# list files
fl() {
  find "${1:-.}" -type f -printf '%T@ %M %n %u %g %14s %Tb %Td %TH:%TM %p\n' | sort -rn | head -${2:-5} | cut -f2- -d" "
}

# nextflow clean
nextflow_clean() {
  nextflow clean -f
  rm -r results.* .nextflow* work
}

weather() {
   curl http://v2.wttr.in
}

# this part is set only if run from an interactive shell
if [[ $- == *i* ]]; then
  
  # color prompt
  export PS1="\[$(tput setaf 7)\]\u@\[$(tput bold)\]\[$(tput setaf 1)\]\h\[$(tput setaf 7)\]:\[$(tput setaf 2)\]\w\[$(tput sgr0)\]\[$(tput setaf 7)\]> \[$(tput sgr0)\]"
  
  # trim prompt length
  export PROMPT_DIRTRIM=3

  # screen prompt
  export PROMPT_COMMAND="prompt_cmd"
  
  # set options for less
  export LESS='--quit-if-one-screen --ignore-case --status-column --LONG-PROMPT --RAW-CONTROL-CHARS --HILITE-UNREAD --tabs=4 --no-init --window=-4'
  export LESSCOLORIZER='pygmentize'

  ### history is saved per-(screen)-session
  # Convert /dev/nnn/X or /dev/nnnX to "nnnX"
  HISTSUFFIX=`tty | sed 's/\///g;s/^dev//g'`
  
  # History file is now .bash_history_pts0
  export HISTFILE="~/.histories/.bash_history_${HISTSUFFIX}_${STY}_$(timestamp)"
  export HISTTIMEFORMAT="%y-%m-%d %H:%M:%S "
  export HISTCONTROL=ignoredups:ignorespace
  shopt -s histappend
  export HISTSIZE=1000000
  export HISTFILESIZE=5000000

fi

# ensure X forwarding is setup correctly, even for screen
XAUTH=~/.Xauthority

if [[ ! -e "${XAUTH}" ]]; then
 # create new ~/.Xauthority file
 xauth
fi

if [[ -z "${XAUTHORITY}" ]]; then
 # export env var if not already available.
 export XAUTHORITY="${XAUTH}"
fi
