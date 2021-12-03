#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M fihorvat@gmail.com
#PBS -m n
#PBS -N pbs.09.construct_supertree.RAxML
#PBS -l select=ncpus=6:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
# input 
INPUT_DIR=.
IN_SEQ=($(find $INPUT_DIR -maxdepth 1 \( -name "*.phy" \)))
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.phy}

# script to run
SCRIPT=/common/WORK/fhorvat/programi/RAxML/standard-RAxML/raxmlHPC-PTHREADS-SSE3

# ----------------Commands------------------- #
# start from the FastTree ML tree
${SCRIPT} -T 6 -p 12467 -m BINCAT -s ${FILE} -n ${BASE}.supertree.raxml

### arguments
# -T = Run on N threads
# -f = Select algorithm
	# The -f option is a very powerful option because, in many cases, it allows you to select what 
	# kind of algorithm RAxML shall execute
	# T option = allows to do a more thorough tree search that uses less lazy, i.e., more  
		# exhaustive SPR moves, in stand-alone mode. This algorithm is typically executed in the  
		# very end  of a search done via ­f a.
	# J option = compute SH­like support values on a given tree passed via ­t. 
		# This option will compute sh-like support values as described here  
		# http://www.ncbi.nlm.nih.gov/pubmed/20525638 on a given tree. The input tree is typically 
		# the best-known ML tree found by a RAxML analysis. Keep in mind that for applying the SH-
		# like test the tree needs to be NNI (Nearest Neighbor Interchange) optimal. Thus, RAxML will 
		# initially try to apply NNI moves to further improve the tree and then compute the SH test 
		# for each inner branch of the tree.
# -p = Random number seed 
	# This allows you to reproduce your results and will help me debug the program. 
	# For all options/algorithms  in RAxML  that require some  sort of randomization,  this option  
	# must be specified.
# -m = Model of Binary (Morphological), Nucleotide, Multi­State, or Amino Acid Substitution
	# input from prottest3
# -t = Specify a user starting tree file name in Newick format
# -s = Specify the name of the alignment data file in PHYLIP or FASTA format
# -n = Specifies the name of the output file
