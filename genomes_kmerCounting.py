import kProcessor as kp
from glob import glob
import os
import sys

if len(sys.argv) < 3:
    sys.exit("run: python genomes_kmerCounting.py <genomes_dir> <kDataframes_output_dir>")

genomes_dir = sys.argv[1]
output_directory = sys.argv[2]

kSize = 21
chunk_size = 10000
hashing_mode = 1 # Integer hashing

genomes = glob(f"{genomes_dir}/*")

for genome_fasta in genomes:
    print(f"counting {genome_fasta} kmers ...")
    kf = kp.kDataFrameMQF(kSize)
    kp.countKmersFromFile(kf,{"mode":hashing_mode},genome_fasta, chunk_size)
    print(f"finished and counted {kf.size()} kmers...")
    print(f"saving ...")
    output_prefix = os.path.join(output_directory, "idx_" + os.path.basename(genome_fasta))
    kf.save(output_prefix)
