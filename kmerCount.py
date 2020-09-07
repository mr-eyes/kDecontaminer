import kProcessor as kp
from glob import glob
import os
import sys

if len(sys.argv) < 3:
    sys.exit("run: python genomes_kmerCounting.py <seq.fa> <kDataframes_output_dir>")

genome_fasta = sys.argv[1]
output_directory = sys.argv[2]

kSize = 21
chunk_size = 10000
hashing_mode = 1 # Integer hashing


kf = kp.kDataFrameMQF(kSize)
kp.countKmersFromFile(kf,{"mode":hashing_mode},genome_fasta, chunk_size)
print(f"finished and counted {kf.size()} kmers...")
print(f"saving ...")
output_prefix = os.path.join(output_directory, "idx_" + os.path.basename(genome_fasta))
kf.save(output_prefix)
