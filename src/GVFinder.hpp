#pragma once

#include <string>
#include <vector>

class GVFinder {
protected:
    std::vector<std::string> sequences;
    std::vector<std::vector<std::string> > clusters;
    std::vector<std::string> results;
    int max_cluster_difference = 12;
public:
    GVFinder(std::string data_path);

    int get_data_size() { return (int)sequences.size(); }
    int get_allele_size() { return (sequences.size() ? (int)sequences[0].size() : 0); }
    
    static bool compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y);
    static int get_max_difference(const std::string &X, const std::string &Y);

    void find_alignment();
    void cluster_msa(const std::vector<std::string> &msa);
    void compare_with_known_results(const std::vector<std::string> &known_results, const std::vector<std::string> &results);
    void output_to_file(const std::vector<std::string> &results, std::string path);
    void calculate_results();
    void output(std::string path);
    void solve();
    int get_max_difference(const std::vector<std::string> &cluster, const std::string &sequence);
    bool belongs_to_cluster(const std::vector<std::string> &cluster, const std::string &sequence, const int &max_cluster_difference);
    std::string get_consensus(const std::vector<std::string> &cluster);
    std::vector<std::string> get_known_results(std::string path="");

};