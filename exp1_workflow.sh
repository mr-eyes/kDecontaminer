# Construct combact de bruign graph of genomes

# Define a constant kmer size
KSIZE=25

# cDBG for genome 1
FASTA=GCF_004143615.1_amil_sf_1.1_genomic.fna
GENOME1_OUTPUT=cDBGk${KSIZE}_${FASTA}
bcalm -kmer-size ${KSIZE} -nb-cores 3 -max-memory 6000 -abundance-min 1 -out ${GENOME1_OUTPUT} -in ${FASTA}

# cDBG for genome 2
FASTA=coral_Mcav.genome_assembly.fasta
GENOME2_OUTPUT=cDBGk${KSIZE}_${FASTA}
bcalm -kmer-size ${KSIZE} -nb-cores 3 -max-memory 6000 -abundance-min 1 -out ${GENOME2_OUTPUT} -in ${FASTA}

# Combine unitigs in one fasta and generate names file for kProcessor
COMBINED_PREFIX="exp1"
python unitigsToKprocessorFormat.py ${COMBINED_PREFIX} ${GENOME1_OUTPUT}.unitigs.fa ${GENOME2_OUTPUT}.unitigs.fa


# Indexing all genomes
python indexing.py ${COMBINED_PREFIX}.fa ${COMBINED_PREFIX}.fa.names