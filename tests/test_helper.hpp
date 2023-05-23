#pragma once

#include <vector>
#include <string>
#include "utilities.hpp"

class Test_helper {
private:
    std::vector<std::string> sample29;
    std::vector<std::string> sample30;

public:
    Test_helper();
    void load_samples();
    int get_number_of_matches(const std::vector<std::string> &results, const std::vector<std::string> &known_results, bool print=false);
    Test_result solve(std::string path_to_file, int method, int sample_id, bool print=false, 
                    int max_cluster_difference=12, int max_size_difference=5);
    bool found_sequence(std::string path_to_file, std::string sequence, int method, int max_cluster_difference, int max_size_difference);
};