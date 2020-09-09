#include "string"
#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include "algorithms.hpp"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <vector>
#include <map>
#include <fstream>

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

std::string base_name(std::string const & path)
{
  return path.substr(path.find_last_of("/\\") + 1);
}

int main(int argc, char ** argv){

    string genome_prefix = argv[1];
    eraseSubStr(genome_prefix, ".mqf");
    cerr << "Loading genome KF ... ";
    auto genomeKF = kDataFrame::load(genome_prefix);
    cerr << "[DONE]" << endl;

    vector<string> samples_files;
    map<string, uint32_t> sample_to_commonKmers;


    for(int i = 2; i < argc; i++){
        string sample_name = argv[i];
        eraseSubStr(sample_name, ".mqf");
        samples_files.emplace_back(sample_name);
    }


    for(auto const & sampleFile : samples_files){
        auto kf = kDataFrame::load(sampleFile);
        cerr << "calculating common kmers ("<< sampleFile <<"&"<< genome_prefix <<") ... ";
        auto commonKmersKF = kProcessor::kFrameIntersect({kf, genomeKF});
        cerr << "[DONE]" << endl;
        uint64_t intersection = kf_size(commonKmersKF);
        sample_to_commonKmers[sampleFile] = intersection;
    }
    
    string genome_basename = base_name(genome_prefix);
    cerr << "Wriring results to " << "contamStats_" + genome_basename + ".tsv" << " ..." << endl;
    ofstream fs;
    string output_file_name = "contamStats_" + genome_basename + ".tsv";
    fs.open(output_file_name);
    fs << "ref";

    for(auto const & sample : sample_to_commonKmers){
        fs << '\t' + sample.first;
    }

    fs << '\n';

    fs << genome_prefix;
    for(auto const & sample : sample_to_commonKmers){
        fs << '\t' + sample.second;
    }
    fs << '\n';
    fs.close();


}