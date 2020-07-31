FASTA_FILE=
grep ">" "${FASTA_FILE}" | cut -c2- |  awk -F' ' '{print $0"\t""${FASTA_FILE}"}' > "${FASTA_FILE}".names