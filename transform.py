import sys
intersection_count= list()

genomes_names = set()
samples_names = set()
samples_kmers = dict()
#with open("allSamples.tsv") as stats:
with open(sys.argv[1]) as stats:
    next(stats)
    for line in stats:
        line = line.strip().split('\t')[:4] # make sure they are only 4
        sample_name, genome_name, common_kmers, sample_kmers = tuple(line)
        print(sample_name, genome_name, common_kmers, sample_kmers)
        intersection_count.append((sample_name, genome_name, int(common_kmers), int(sample_kmers)))
        genomes_names.add(genome_name)
        samples_names.add(sample_name)
        samples_kmers[sample_name] = int(sample_kmers)

genomes_names = list(genomes_names)
samples_names = list(samples_names)

intersection_by_genome = dict()

for _genome in genomes_names:
    intersection_by_genome[_genome] = dict()

for item in intersection_count:
    sample_name, genome_name, common_kmers, unused_sample_kmers = item
    intersection_by_genome[genome_name][sample_name] = common_kmers

print(intersection_by_genome)

with open("transformed_" + sys.argv[1], 'w') as OUT:
    header = str()
    header += "ref.\t"
    for sample_name, common_kmers in intersection_by_genome[genomes_names[0]].items():
        sample = sample_name.replace(".mqf", '')
        header += sample + '\t'
    OUT.write(header[:-1] + '\n')

    print(sample_kmers)

    for genome_name, sample_dict in intersection_by_genome.items():
        row = f"{genome_name}\t"
        for sample_name, common_kmers in sample_dict.items():
            containment = 100 * (common_kmers / samples_kmers[sample_name])
            row += '%' + str(containment) + '\t'

        OUT.write(row[:-1] + '\n')

print("samples kmers:")
for sample_name, kmers in samples_kmers.items():
    print(f"\tsample({sample_name}): {kmers} kmers")