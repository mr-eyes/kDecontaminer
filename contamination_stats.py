import kProcessor as kp
from glob import glob
import os
import sys
import multiprocessing as MP
import pickle
import gc

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

_sample_counter = 1

for sample in samples_kfs:
    sample = sample.replace(".mqf", '')
    print(f"Getting sample {sample} ({_sample_counter}) kmers")
    tmp_kf = kp.kDataFrame.load(sample)
    sample_kmers[os.path.basename(sample)] = tmp_kf.size()
    samples_names.append(os.path.basename(sample))

    for genome in genome_kfs:
        genome = genome.replace(".mqf", '')
        genomes_names.append(os.path.basename(genome))
        job_pairs.append((sample, genome))

manager = MP.Manager()

intersection_count = manager.list()

preloaded_kfs = manager.dict()

for genome in genome_kfs:
    genome = genome.replace(".mqf", '')
    print(f"loading {genome}")
    preloaded_kfs[genome] = kp.kDataFrame.load(genome)

for sample in samples_kfs:
    sample = sample.replace(".mqf", '')
    print(f"loading {sample}")
    preloaded_kfs[sample] = kp.kDataFrame.load(sample)


def get_intersection(pair):
    global intersection_count
    global preloaded_kfs

    print(f"processing ({pair})")

    sample_file = pair[0]
    genome_file = pair[1]

    intersection_kf = kp.kFrameIntersect([preloaded_kfs[sample_file], preloaded_kfs[genome_file]])

    common_kmers = intersection_kf.size()
    sample_kmers = preloaded_kfs[sample_file].size()

    sample_name = os.path.basename(pair[0])
    genome_name = os.path.basename(pair[1])

    intersection_count.append((sample_name, genome_name, common_kmers, sample_kmers))
    print("---------------------------------")
    print(pair)
    print(sample_name, genome_name, common_kmers, sample_kmers)
    print("---------------------------------")

    del intersection_kf
    gc.collect()

print(f"Processing started ...")

with MP.Pool(threads) as pool:
    pool.map(get_intersection, job_pairs)
    pool.close()
    pool.join()

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
