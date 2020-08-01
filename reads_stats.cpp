#include "colored_kDataFrame.hpp"
#include <string>
#include <iostream>
#include <vector>
#include "kseqReader.hpp"
#include <sys/stat.h>
#include <map>
#include "readClassifier.hpp"
#include <cassert>


using namespace std;

inline bool file_exists(const std::string &name) {
    struct stat buffer
            {
            };
    return (stat(name.c_str(), &buffer) == 0);
}

inline string time_diff(std::chrono::high_resolution_clock::time_point &t1) {
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long hr = milli / 3600000;
    milli = milli - 3600000 * hr;
    long min = milli / 60000;
    milli = milli - 60000 * min;
    long sec = milli / 1000;
    milli = milli - 1000 * sec;
    string timeDiff;
    timeDiff.append(to_string(min));
    timeDiff.append(":");
    timeDiff.append(to_string(sec));
    timeDiff.append(":");
    timeDiff.append(to_string(milli));

    return timeDiff;
}


class stats {


    stats(int no_genomes) {

    }

};


int main(int argc, char **argv) {

    if (argc != 4) {
        cerr << "run ./peReadsStats <R1> <R2> <index_prefix>" << endl;
        exit(1);
    }

    string PE_1_reads_file = argv[1];
    string PE_2_reads_file = argv[2];
    string index_prefix = argv[3];

    if (!file_exists(PE_1_reads_file)) {
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(PE_2_reads_file)) {
        throw std::runtime_error("Could not open R2 file");
    }

    if (!file_exists(index_prefix + ".extra")) {
        throw std::runtime_error("Could not open kProcessor index file");
    }

    int batchSize = 5000;
    int kSize = 25;
    int no_of_sequences;

    cerr << "counting number of reads ..." << endl;
    int count = 0;
    string line;
    ifstream file(PE_1_reads_file);
    while (getline(file, line)) count++;

    if (PE_1_reads_file.find("fastq") != std::string::npos || PE_1_reads_file.find("fq") != std::string::npos) {
        no_of_sequences = count / 4;
    }else{
        no_of_sequences = count / 2;
    }

    int no_chunks = ceil((double) no_of_sequences / (double) batchSize);
    cerr << "processing " << no_of_sequences << " reads in " << no_chunks << " chunks ..." << endl;

    // kProcessor Index Loading
    std::cerr << "Loading kProcessor index ..." << std::endl;
    colored_kDataFrame *ckf = colored_kDataFrame::load(index_prefix);
    kDataFrame *kf = ckf->getkDataFrame();

    flat_hash_map<uint64_t, int> singleColors;
    vector<int> vec_singleColors;
    for (auto const &color: ckf->namesMap) {
        singleColors[color.first] = 1;
        vec_singleColors.push_back(color.first);
    }

    assert(kSize == (int) kf->getkSize());
    assert(kf->size() > 100);
    std::cerr << "kProcessor index loaded successfully ..." << std::endl;

    // Initializations

    int current_chunk = 0;

    auto *PEReader = new kseqReader(PE_1_reads_file, PE_2_reads_file, batchSize);

    auto *tmpReader = new kseqReader(PE_1_reads_file, PE_2_reads_file, 1);
    unsigned long read_length = tmpReader->next_chunk()->back().R1_seq.size();
    delete tmpReader;


    // Output TSV file
    ofstream output("stats_" + index_prefix.substr(index_prefix.find_last_of("/\\") + 1) + ".tsv");
    for (auto const &color: vec_singleColors) {
        output << "uniq_genome_" << color << '\t';
    }
    output << "ambiguous" << '\t' << "total_aligned" << '\t' << "unmapped" << '\n';


    while (!PEReader->end()) {
        std::chrono::high_resolution_clock::time_point _currentTime = std::chrono::high_resolution_clock::now();
        cerr << "Processing chunk " << ++current_chunk << "/" << no_chunks << "... ";
        vector<tuple<int, string, string, string>> sqlite_chunk; // Buffer for holding Sqlite rows

        for (auto const &PE : *PEReader->next_chunk()) {
            vector<string> kmers;


            flat_hash_map<int, int> uniqueCount;
            for (auto const &c : vec_singleColors) {
                uniqueCount[c] = 0;
            }
            int ambiguous = 0 , unmapped = 0;

            for (unsigned long i = 0; i < read_length - kSize + 1; i++) {
                kmers.emplace_back(PE.R1_seq.substr(i, kSize));
                kmers.emplace_back(PE.R2_seq.substr(i, kSize));
            }

            for (auto const &kmer : kmers) {

                // Exclude any kmer with Ns
                if (kmer.find('N') != std::string::npos) {
                    // cout << kmer << endl;
                    unmapped++;
                    continue;

                } else {
                    // Get the kmer color
                    uint64_t color = kf->getCount(kmer);

                    // Check if the alignment is unique to one genome
                    if (singleColors[color]) {
                        uniqueCount[color]++;
                    } else if (color) {
                        ambiguous++;
                    } else {
                        unmapped++;
                    }
                }

            }

            int total_aligned = 0;

            for (auto const &color: vec_singleColors) {
                output << uniqueCount[color] << '\t';
                total_aligned += uniqueCount[color];
            }
            output << ambiguous << '\t' << total_aligned + ambiguous << '\t' << unmapped << '\n';
        }

        cerr << "Done in " << time_diff(_currentTime) << endl;
    }


    return 0;
}
