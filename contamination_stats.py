import kProcessor as kp
from glob import glob
import os
import sys
import multiprocessing as MP
import pickle

if len(sys.argv) < 5:
    sys.exit("run: python genomes_kmerCounting.py <genomes_KFs_dir> <reads_cDBGs_dir> <output_file> <threads>")

genomes_dir = sys.argv[1]
reads_dir = sys.argv[2]
output_file = sys.argv[3] + ".tsv"
threads = int(sys.argv[4])

kSize = 21
chunk_size = 10000
hashing_mode = 1  # Integer hashing

samples_kfs = glob(reads_dir + "/*.mqf")
genome_kfs = glob(genomes_dir + "/*.mqf")

job_pairs = list()

genomes_names = list()
samples_names = list()
sample_kmers = dict()

for sample in samples_kfs:
    sample = sample.replace(".mqf", '')
    tmp_kf = kp.kDataFrame.load(sample)
    sample_kmers[os.path.basename(sample)] = tmp_kf.size()
    samples_names.append(os.path.basename(sample))

    for genome in genome_kfs:
        genome = genome.replace(".mqf", '')
        genomes_names.append(os.path.basename(genome))
        job_pairs.append((sample, genome))

manager = MP.Manager()

intersection_count = manager.list()

def get_intersection(pair):
    global intersection_count
    print(f"processing ({pair})")
    sample_kf = kp.kDataFrame.load(pair[0])
    genome_kf = kp.kDataFrame.load(pair[1])

    intersection_kf = kp.kFrameIntersect([sample_kf, genome_kf])

    common_kmers = intersection_kf.size()
    sample_kmers = sample_kf.size()

    sample_name = os.path.basename(pair[0])
    genome_name = os.path.basename(pair[1])

    intersection_count.append((sample_name, genome_name, common_kmers, sample_kmers))


print(f"Processing started ...")
with MP.Pool(threads) as pool:
    pool.map(get_intersection, job_pairs)

with open(f"{output_file}.pickle", "wb") as fp:
    pickle.dump(intersection_count, fp)
print("dumped the result in pickle ...")

intersection_by_genome = dict()

for _genome in genomes_names:
    intersection_by_genome[_genome] = dict()

for item in intersection_count:
    sample_name, genome_name, common_kmers, unused_sample_kmers = item
    intersection_by_genome[genome_name][sample_name] = common_kmers

with open(output_file, 'w') as OUT:
    header = str()
    header += "ref.\t"
    for sample_name, common_kmers in intersection_by_genome[genomes_names[0]].items():
        sample = sample_name.replace(".mqf", '')
        header += sample + '\t'
    OUT.write(header[:-1] + '\n')

    for genome_name, sample_dict in intersection_by_genome.items():
        row = f"{genome_name}\t"
        for sample_name, common_kmers in sample_dict.items():
            containment = 100 * (common_kmers / sample_kmers[sample_name])
            row += '%' + str(containment) + '\t'

        OUT.write(row[:-1] + '\n')

print("samples kmers:")
for sample_name, kmers in sample_kmers.items():
    print(f"\tsample({sample_name}): {kmers} kmers")
