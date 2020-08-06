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

inline bool file_exists(const std::string &name)
{
    struct stat buffer
            {
            };
    return (stat(name.c_str(), &buffer) == 0);
}

int main(int argc, char **argv) {

    if (argc != 2) {
        cerr << "run ./analyzeIdx <index_prefix>" << endl;
        exit(1);
    }

    string index_prefix = argv[1];

    if (!file_exists(index_prefix + ".extra")) {
        throw std::runtime_error("Could not open kProcessor index file");
    }

    colored_kDataFrame *ckf = colored_kDataFrame::load(index_prefix);
    kDataFrame *kf = ckf->getkDataFrame();


    flat_hash_map<uint64_t, int> singleColors;
    vector<int> vec_singleColors;
    for (auto const &color: ckf->namesMap) {
        singleColors[color.first] = 1;
        vec_singleColors.push_back(color.first);
    }

    cout << "Single colors:" << endl;


    cout << "Index loaded successfully .." << endl;

    auto kIt = kf->begin();

    flat_hash_map<int, uint64_t> stats;


    while(kIt != kf->end()){
        stats[kIt.getCount()]++;
        kIt++;
    }

    for(auto const & color: stats){
        cout << "color(" << color.first <<"): " << color.second << endl;
    }

}