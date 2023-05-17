#pragma once

#include <string>
#include <vector>

class GVFinder {
protected:
    std::vector<std::string> sequences;
    std::vector<std::vector<std::string> > clusters;
    int max_cluster_difference = 12;
public:
    GVFinder(std::string data_path);

    int get_data_size() { return (int)sequences.size(); }
    int get_allele_size() { return (sequences.size() ? (int)sequences[0].size() : 0); }
    
    static bool compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y);

    void find_alignment();
    void cluster_msa(const std::vector<std::string> &msa);
    int get_max_difference(const std::vector<std::string> &cluster, const std::string &sequence);
};