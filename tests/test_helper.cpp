#include "test_helper.hpp"
#include "utilities.hpp"
#include <fstream>
#include <string>
#include "../src/GVFinder.hpp"
#include <assert.h>
#include <iostream>

void Test_helper::load_samples() {
    std::ifstream stream;
    stream.open("data/alela_jelen_modified.fasta");
    std::string info;
    std::string solution;
    while(std::getline(stream, info)){
        if(std::getline(stream, solution)){
            if(info.rfind(">J30B", 0) == 0){
                sample30.push_back(solution);
            } else if(info.rfind(">J29B") == 0){
                sample29.push_back(solution);
            }
        }
    }

    for(auto &sample : sample29)
        while((int)sample.size() > 249)
            sample.pop_back();

    for(auto &sample : sample30)
        while((int)sample.size() > 249)
            sample.pop_back();
    stream.close();
}

Test_helper::Test_helper() {
    load_samples();
}

int Test_helper::get_number_of_matches(const std::vector<std::string> &results, const std::vector<std::string> &known_results, bool print){
    if(print) {
        std::cout << "Results:\n";
        for(const auto &it : results)
            std::cout << it << "\n";
        std::cout << "Known results:\n";
        for(const auto &it : known_results)
            std::cout << it << "\n";
        std::cout << "\n";
    }
    
    int number_of_matches = 0;
    for(const auto &known_result : known_results){
        for(const auto &result : results){
            if(result.find(known_result) != std::string::npos){
                number_of_matches ++;
                break;
            }
        }
    }
    return number_of_matches;
}

Test_result Test_helper::solve(std::string path_to_file, int method, int sample_id, bool print,
                                int max_cluster_difference, int max_size_difference){
    GVFinder gvfinder = GVFinder(path_to_file, method);
    gvfinder.set_max_cluster_difference(max_cluster_difference);
    gvfinder.set_max_size_difference(max_size_difference);

    gvfinder.solve();
    const std::vector<std::string> results = gvfinder.get_results();
    if(sample_id == 29)
        return {get_number_of_matches(results, sample29, print), (int)results.size()};
    else if(sample_id == 30)
        return {get_number_of_matches(results, sample30, print), (int)results.size()};
    else
        assert(false);
    return {-1, -1};
}