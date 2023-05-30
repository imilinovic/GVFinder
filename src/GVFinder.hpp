#pragma once

#include <string>
#include <vector>

class GVFinder {
protected:
    std::vector<std::string> sequences;
    std::vector<std::vector<std::string> > clusters;
    std::vector<std::vector<std::string> > clusters2;
    std::vector<std::string> results;
    int max_cluster_difference;
    int max_size_difference;
    int method;

public:
    GVFinder(std::string data_path, int method, int max_size_difference, int max_cluster_difference);

    int get_data_size() { return (int)sequences.size(); }
    int get_allele_size() { return (sequences.size() ? (int)sequences[0].size() : 0); }
    void set_max_cluster_difference(int max_cluster_difference) { this -> max_cluster_difference = max_cluster_difference; }
    void set_max_size_difference(int max_size_difference) { this -> max_size_difference = max_size_difference; }
    const std::vector<std::string> &get_results() { return results; }

    static bool compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y);
    static int get_max_difference(const std::string &X, const std::string &Y);

    void find_alignment();
    void cluster_msa_1(const std::vector<std::string> &msa);
    void cluster_msa_2(const std::vector<std::string> &msa);
    void cluster_msa_3(const std::vector<std::string> &msa);
    void compare_with_known_results(const std::vector<std::string> &known_results, const std::vector<std::string> &results);
    void output_to_file(std::string filename, const std::vector<std::string> &results, std::string path);
    void calculate_results();
    void output(std::string filename, std::string path="");
    void check(std::string ssequence, std::string name);
    void solve();
    int get_max_difference(const std::vector<std::string> &cluster, const std::string &sequence);
    bool belongs_to_cluster(const std::vector<std::string> &cluster, const std::string &sequence, const int &max_cluster_difference);
    double calculate_average_distance(const std::vector<std::string> &cluster, const std::string &sequence, const int &max_cluster_difference);
    std::string get_consensus(const std::vector<std::string> &cluster);
    std::vector<std::string> get_known_results(std::string path="");
};