#include "string"
#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <vector>

using namespace std;

void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

uint64_t kf_size(kDataFrame * KF){

    uint64_t size = 0;
    auto it = KF->begin();
    while(it != KF->end()){
        size++;
        it++;
    }
    return size;
}

int main(int argc, char ** argv){

    string genome_prefix = argv[1];
    eraseSubStr(genome_prefix, ".mqf");

    auto genomeKF = kDataFrame::load(genome_prefix);

    cout << "size(" << genome_prefix << ") = " << kf_size(genomeKF) << endl;

}