import os
import sys

import kProcessor as kp

fasta_file = str()
namesFile = str()

if len(sys.argv) < 4:
    sys.exit("run: python genes_indexing.py <fasta> <namesFiles> <kSize>")

else:
    fasta_file = sys.argv[1]
    namesFile = sys.argv[2]
    kSize = int(sys.argv[3])


idx_suffix = os.path.basename(fasta_file).replace(".fa.names", "")
print(f"Indexing {idx_suffix} ...", file=sys.stderr)
kmers_mode = 1
hashing_mode = 1
chunk_size = 10000
#kf_PHMAP = kp.kDataFramePHMAP(kSize, hashing_mode)
kf_PHMAP = kp.kDataFrameMQF(kSize, 29, hashing_mode)
ckf = kp.index(kf_PHMAP, {"mode": kmers_mode}, fasta_file, chunk_size, namesFile)
ckf.save(f"idx_{idx_suffix}")
