#include <cstdint>
#include "colored_kDataFrame.hpp"
#include <parallel_hashmap/phmap.h>

class Statistics
{
public:
    uint64_t n_unmatched;
    uint64_t n_same;
    uint64_t n_amb_same;
    uint64_t n_fusion;
    uint64_t n_ambig_fusion;
    uint64_t n_mutli_fusion;
    uint64_t n_paired_fusion;
    uint64_t n_paired_unique;
    uint64_t n_fragments;
    uint64_t n_paired_ambig_unique;
    uint64_t n_paired_unmatched;

    Statistics()
    {
        n_fragments = 0;
        n_unmatched = 0;
        n_same = 0;
        n_amb_same = 0;
        n_fusion = 0;
        n_ambig_fusion = 0;
        n_mutli_fusion = 0;
        n_paired_fusion = 0;
        n_paired_unique = 0;
        n_paired_ambig_unique = 0;
        n_paired_unmatched = 0;
    }

    void print() const
    {
        cout << "Summary report" << endl;
        cout << "===============" << endl
             << endl;
        cout << "Total fragments: " << n_fragments << endl
             << endl;

        cout << "Reads stats:" << endl;

        cout << "\tunique:             " << n_same << endl;
        cout << "\tambiguous:          " << n_amb_same << endl;
        cout << "\tchimeric:           " << n_fusion << endl;
        cout << "\tunmatched:          " << n_unmatched << endl;

        cout << endl;

        cout << "Fragments stats:" << endl;

        cout << "\tunique:             " << n_paired_unique << endl;
        cout << "\tambiguous:          " << n_paired_ambig_unique << endl;
        cout << "\tchimeric:           " << n_paired_fusion << endl;
        cout << "\tunmatched:          " << n_paired_unmatched << endl;
    }
};

class ReadsClassifier
{

private:
    colored_kDataFrame *DB;
    int kSize;
    unordered_map<int, string> namesMap;

    vector<vector<uint32_t>> families0, families1, families;
    vector<uint32_t> shared_kmers0, gaps0, shared_kmers1, gaps1, shared_kmers, gaps;

public:
    Statistics stats;

    explicit ReadsClassifier(colored_kDataFrame *ckf)
    {
        this->DB = ckf;
        this->kSize = ckf->getkDataFrame()->ksize();
        namesMap = ckf->names_map();
    }

    static bool isSameRef(string &flag)
    {
        return flag == "unique" || flag == "ambiguous";
    }

    template <class T>
    static string setsToString(std::vector<T> &v)
    {
        if (v.size() == 0)
        {
            return "{}";
        }
        string res = "{" + to_string(v[0]);
        for (int i = 1; i < v.size(); i++)
        {
            res += " ," + std::to_string(v[i]);
        }
        res += "}";
        return res;
    }

    template <class T>
    static string familiesToString(vector<vector<T>> &v)
    {
        if (v.size() == 0)
        {
            return "[]";
        }
        string res = "[" + setsToString(v[0]);
        for (int i = 1; i < v.size(); i++)
        {
            res += " ," + setsToString(v[i]);
        }
        res += "]";
        return res;
    }

    tuple<string, string, string, vector<uint32_t>> classifyFragment(const string &R1_seq, const string &R2_seq)
    {

        families0.clear();
        gaps0.clear();
        shared_kmers0.clear();
        families1.clear();
        gaps1.clear();
        shared_kmers1.clear();

        vector<uint32_t> final_result_vector;

        string flag_R1 = readFusion(R1_seq, families0, shared_kmers0, gaps0);
        string flag_R2 = readFusion(R2_seq, families1, shared_kmers1, gaps1);

        // cout << flag_R1 << '\t' << familiesToString(families0) << endl;
        // cout << flag_R2 << '\t' << familiesToString(families1) << endl;
        // return make_tuple("","","",final_result_vector);

        string fragmentFlag{};

        /*
         * unique, ambiguous, unmatched
         * multi_fusion, ambig_fusion, clear_fusion
         * */

        vector<uint32_t> partitions_ids;

        if (flag_R1 == "fusion" || flag_R2 == "fusion")
        {
            // One or both are chimeric
            fragmentFlag = "chimeric";
            stats.n_paired_fusion++;
        }
        else if (isSameRef(flag_R1) && isSameRef(flag_R2))
        {
            // Both are ambiguous or unique but may be no intersection
            partitions_ids.clear();
            auto it = set_intersection(families0[0].begin(), families0[0].end(), families1[0].begin(), families1[0].end(), back_inserter(partitions_ids));

            if (partitions_ids.empty())
            {
                // If no intersection, then they are chimeric
                stats.n_paired_fusion++;
                fragmentFlag = "chimeric";
            }
            else
            {
                // Check if they are unique or ambiguous
                if (flag_R1 == "unique" && flag_R2 == "unique")
                {
                    fragmentFlag = "unique";
                    stats.n_paired_unique++;
                }
                else
                {
                    if (partitions_ids.size() == 1)
                    {
                        // each read is marked as ambiguous but the intersection is unique
                        fragmentFlag = "unique";
                        stats.n_paired_unique++;
                    }
                    else // Intersection is > 1 then it's an obvious ambiguous
                    {
                        fragmentFlag = "ambiguous";
                        stats.n_paired_ambig_unique++;
                    }
                }
            }
        }
        else
        {

            if (isSameRef(flag_R1) || isSameRef(flag_R2))
            {
                // Reads are ambiguous or unique

                if (flag_R1 == "unmatched" || flag_R2 == "unmatched")
                {
                    partitions_ids.clear();

                    if (flag_R1 == "unmatched")
                    {
                        families0[0].clear();
                        partitions_ids = families1[0];
                    }
                    else
                    {
                        families1[0].clear();
                        partitions_ids = families0[0];
                    }
                    // One of the reads must be unmatched
                    // Change partition_ids to union instead of intersect
                    // partitions_ids.erase(std::remove(partitions_ids.begin(), partitions_ids.end(), 0), partitions_ids.end());

                    auto end = partitions_ids.end();
                    for (auto it = partitions_ids.begin(); it != end; ++it)
                        end = std::remove(it + 1, end, *it);
                    partitions_ids.erase(end, partitions_ids.end());

                    if (flag_R1 == "unique" || flag_R2 == "unique")
                    {

                        fragmentFlag = "unique";
                        stats.n_paired_unique++;
                    }

                    else if (flag_R1 == "ambiguous" || flag_R2 == "ambiguous")
                    {
                        // One of them must be ambiguous

                        fragmentFlag = "ambiguous";
                        stats.n_paired_ambig_unique++;
                    }
                }
            }
            else
            {
                fragmentFlag = "unmatched";
                stats.n_paired_unmatched++;
            }
        }

        return make_tuple(flag_R1, flag_R2, fragmentFlag, partitions_ids);
    }

    string readFusion(const string &read, vector<vector<uint32_t>> &families,
                      vector<uint32_t> &shared_kmers, vector<uint32_t> &gaps)
    {

        shared_kmers.clear();
        families.clear();
        gaps.clear();
        string flag;
        vector<uint32_t> lf_ids;
        vector<uint32_t> rt_ids;

        if (read.size() < kSize)
        {
            stats.n_unmatched++;
            return "unmatched";
        }
        //# find a matching k-mer at the beginning of the read

        lf_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(0, kSize)));

        int idx = 1;

        while (idx < read.size() - kSize + 1 && lf_ids.empty())
        {
            lf_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(idx, kSize)));
            idx++;
        }
        if (lf_ids.empty())
        {
            //        print("no single match");
            stats.n_unmatched++;
            flag = "unmatched";
        }
        else if (idx == read.size() - kSize + 1)
        {
            //        print("same, only last kmer matched");
            families.push_back(lf_ids);

            if (lf_ids.size() == 1)
            {
                stats.n_same += 1;
                flag = "unique";
            }
            else
            {
                stats.n_amb_same += 1;
                flag = "ambiguous";
            }
        }
        else
        {
            // # len(lf_ids) > 0 & idx < len(hashvals)
            //# find a matching k-mer at the end of the read
            vector<uint32_t> rt_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(read.size() - kSize, kSize)));

            int idy = read.size() - 1;
            while (idy - kSize >= idx - 1 && rt_ids.empty())
            {
                rt_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(idy - kSize, kSize)));
                idy--;
            }

            if (rt_ids.empty())
            {
                //            print("same, only one non-last kmer matched");
                families.push_back(lf_ids);
                if (lf_ids.size() == 1)
                {
                    stats.n_same += 1;
                    flag = "unique";
                }
                else
                {
                    stats.n_amb_same += 1;
                    flag = "ambiguous";
                }
            }
            else
            {
                vector<uint32_t> partitions_ids;
                partitions_ids.clear();
                auto it = set_intersection(lf_ids.begin(), lf_ids.end(), rt_ids.begin(), rt_ids.end(), back_inserter(partitions_ids));
                if (!partitions_ids.empty())
                {
                    families.push_back(partitions_ids);
                    if (partitions_ids.size() == 1)
                    {
                        stats.n_same += 1;
                        flag = "unique";
                    }
                    else
                    {
                        stats.n_amb_same += 1;
                        flag = "ambiguous";
                    }
                }
                else
                {
                    //# fusion to be resolved
                    flag = "fusion";
                    stats.n_fusion += 1;
                }
            }
        }
        return flag;
    }
};