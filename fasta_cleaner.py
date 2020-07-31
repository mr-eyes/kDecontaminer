import os
import sys

if len(sys.argv) < 2:
    sys.exit("run: python indexing.py <fasta_file>")
else:
    fasta_file = sys.argv[1]

if not os.path.exists(fasta_file):
    sys.exit(f"{fasta_file} file does not exist.")

clean_fasta = "clean_" + os.path.basename(fasta_file)

with open(fasta_file) as FASTA, open(clean_fasta, 'w') as CLEAN_FASTA:
    for line in FASTA:
        if ">" in line:
            line = line.replace('\t', ' ')
            CLEAN_FASTA.write(line)
        else:
            CLEAN_FASTA.write(line)