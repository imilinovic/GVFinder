#include "test_helper.hpp"
#include "utilities.hpp"
#include <fstream>
#include <string>
#include "../src/GVFinder.hpp"
#include <assert.h>

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

int Test_helper::get_number_of_matches(const std::vector<std::string> &results, const std::vector<std::string> &known_results){
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

Test_result Test_helper::solve(std::string path_to_file, int method, int sample_id){
    GVFinder gvfinder = GVFinder(path_to_file, method);
    gvfinder.solve();
    const std::vector<std::string> results = gvfinder.get_results();
    if(sample_id == 29)
        return {get_number_of_matches(results, sample29), (int)results.size()};
    else if(sample_id == 30)
        return {get_number_of_matches(results, sample30), (int)results.size()};
    else
        assert(false);
}