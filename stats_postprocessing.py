import os
import sys

if len(sys.argv) < 2:
    sys.exit("run: python stats_postprocessing.py <stats.tsv>")

total = dict() # {genome_id : indexes}

count = dict()

with open(sys.argv[1]) as statsReader:
    header = next(statsReader).strip().split()
    items_no = len(header)
    all_genomes = list()
    for item in range(len(header)):
        if header[item].startswith("uniq"):
            genome_id = int(header[item].split("_")[-1])
            total[genome_id] = set([item])
            count[genome_id] = 0
            all_genomes.append(genome_id)
        elif header[item].startswith("ambig"):
            genome_ids = list(map(int, header[item].replace("ambig(", "").replace(")","").split("-")))
            for gID in genome_ids:
                total[gID].add(item)


    for line in statsReader:
        line = list(map(int, line.strip().split()))

        for genome in all_genomes:
            for idx in total[genome]:
                count[genome] += line[idx]            


total_count = sum(count.values())


for genome in all_genomes:
    print(f"Genome({genome}): {count[genome]} matched kmers | ~{(100*count[genome]/total_count):.2f}%")

