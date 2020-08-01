#include "kseqReader.hpp"


kseqReader::kseqReader(string R1File, string R2File, int chunk_size) {
    this->chunk_size = chunk_size;
    this->R1_fastx = R1File;
    this->R2_fastx = R2File;

    this->fp_r1 = gzopen(R1File.c_str(), "r");
    this->fp_r2 = gzopen(R2File.c_str(), "r");
    this->seq_1 = kseq_init(this->fp_r1);
    this->seq_2 = kseq_init(this->fp_r2);
}


std::vector<peRead> *kseqReader::next_chunk() {
    this->chunk_reads.clear();

    for (int i = 0; i < this->chunk_size && ((kseq_read(this->seq_1)) >= 0 && (kseq_read(this->seq_2)) >= 0); i++) {
        this->chunk_reads.emplace_back(peRead{this->seq_1->seq.s, this->seq_2->seq.s, this->seq_1->name.s, this->seq_2->name.s});
    }

    if ((int) this->chunk_reads.size() != this->chunk_size) {
        this->END = true;
    }

    return &this->chunk_reads;

}