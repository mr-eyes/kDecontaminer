import kProcessor as kp
from glob import glob
import os
import sys

if len(sys.argv) < 3:
    sys.exit("run: python genomes_kmerCounting.py <genomes_KFs_dir> <reads_cDBGs_dir>")

genomes_dir = sys.argv[1]
reads_dir = sys.argv[2]

kSize = 21
chunk_size = 10000
hashing_mode = 1 # Integer hashing

reads_files = glob("f{reads_dir}/*")

samples_kfs = dict()
samples_total_kmer_count = dict()


for readsFile in reads_files:
    file_name = os.path.basename(readsFile)
    samples_kfs[file_name] = kp.kDataFrameMQF(kSize)
    kp.countKmersFromFile(samples_kfs[file_name] ,{"mode":hashing_mode},readsFile, chunk_size)
    samples_total_kmer_count[file_name] = samples_kfs[file_name].size()


genomes_kfs = glob(f"{genomes_dir}/idx_*.mqf")

intersection_count = dict()

for genome_kf in genomes_kfs:
    kf_prefix = genome_kf.replace(".mqf","")
    
    print(f"Loading genome: {kf_prefix}")
    genome_kf = kp.kDataFrame.load(kf_prefix)
    

    for sample_name, readsKF in samples_kfs.items():
        print(f"Processing read({sample_name} with genome({kf_prefix}))")
        intersection_kf = kp.kFrameIntersect([readsKF, genome_kf])
        intersection_count[tuple([sample_name, kf_prefix])] = intersection_kf.size()


for record, count in intersection_count.items():
    sample_name = record[0]
    genome_name = record[1]
    sample_kmers_count = samples_total_kmer_count[sample_name]
    common_kmers = count
    percentage = 100 * (count / sample_kmers_count)
    print(f"Sample ({sample_name}) has ~{percentage:.2f}% matched kmers on Genome {genome_name}")