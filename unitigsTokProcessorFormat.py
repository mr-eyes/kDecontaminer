from tqdm import tqdm
import subprocess
import sys
import os


if len(sys.argv) < 3:
    sys.exit("run: python unitigsTokProcessorFormat.py <output_prefix> <fasta1> <fasta2> ...")

output_prefix = sys.argv[1]

with open(output_prefix + ".fa" , 'w') as FASTA, open(output_prefix + ".fa.names", 'w') as NAMES:
    for unitigsFile in sys.argv[2:]:
        print(f"Processing {unitigsFile} ...")
        no_seqs = int(subprocess.getoutput('wc -l ' + unitigsFile).split()[0]) // 2
        genome_name = os.path.basename(unitigsFile)
        genome_name = genome_name.replace(".unitigs.fa","")
        genome_name = "".join(genome_name.split('.')[:-1])
        unitigID = 1
        with open(unitigsFile, 'r') as unitigsReader:
            for line in tqdm(unitigsReader, total=no_seqs):
                header = f">{genome_name}_{unitigID}\n"
                seq = next(unitigsReader)
                FASTA.write(header+seq)
                NAMES.write(f"{genome_name}_{unitigID}\t{genome_name}\n")   
                unitigID += 1             