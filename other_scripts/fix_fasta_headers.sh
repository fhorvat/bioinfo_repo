# replaces ":" and "-" with "_"
sed -E -e '/^>/s/-|:/_/g' ${FASTA_FILE} > ${FASTA_FILE_FIXED}
