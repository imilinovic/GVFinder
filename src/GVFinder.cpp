#include "GVFinder.hpp"
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include "spoa/spoa.hpp"


GVFinder::GVFinder(std::string data_path, int method) {
    std::vector<std::string> data;
    std::ifstream stream;
    std::unordered_map<int, int> length_cnt;
    this->method = method;
    



    std::string method_name;
    switch(method){
        case 1:
            method_name = "big punishment";
            break;
        case 2:
            method_name = "rough";
            break;
        case 3:
            method_name = "closest cluster";
            break;
        default:
            throw std::invalid_argument("Invalid method argument");
            break;
    }

    std::cout << "Method: " << method_name << "\n";

    stream.open(data_path);
    std::string line;
    while(std::getline(stream, line)){
        if(std::getline(stream, line)){
            data.push_back(line);
            int len = (int)line.size();
            if(length_cnt.find(len) != length_cnt.end())
                length_cnt[len] ++;
            else
                length_cnt[len] = 1;
        }
        std::getline(stream, line);
        std::getline(stream, line);
    }
    stream.close();

    int most_common_len = -1, max_len_cnt = 0;
    for(auto len : length_cnt){
        if(len.second > max_len_cnt){
            max_len_cnt = len.second;
            most_common_len = len.first;
        }
    }

    if(method == 1){
        for(auto allele : data){
            if((int)allele.size() == most_common_len){
                sequences.push_back(allele);
            }
        }
    }else{
        for(auto allele : data){    
            if((int)allele.size() <= most_common_len + 5 && (int)allele.size() >= most_common_len - 5){
                sequences.push_back(allele);
            }
        }
    }
}

void GVFinder::find_alignment() {
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 0, -1, method == 1 ? -100 : -1);  // linear gaps match mismatch gap
    spoa::Graph graph{};

    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    auto consensus = graph.GenerateConsensus();

    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
            << consensus << std::endl;

    auto msa = graph.GenerateMultipleSequenceAlignment();
    std::cerr << "Velicina MSA: " << (int)msa.size() << "\n";

    switch(method){
        case 1:
            cluster_msa_1(msa);
            break;
        case 2:
            cluster_msa_2(msa);
            break;
        case 3:
            cluster_msa_3(msa);
    }
}

bool GVFinder::compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y){
    return X.size() > Y.size();
}

// int GVFinder::get_max_difference(const std::vector<std::string> &cluster, const std::string &sequence){
//     int max_difference = 0;
//     for(const auto &cluster_sequence : cluster){
//         int difference = 0;
//         for(int i = 0; i < (int)sequence.size(); ++i)
//             if(sequence[i] != cluster_sequence[i])
//                 difference ++;
//         max_difference = std::max(max_difference, difference);
//     }
//     return max_difference;
// }

bool GVFinder::belongs_to_cluster(const std::vector<std::string> &cluster, const std::string &sequence, const int &max_cluster_difference){
    for(const auto &cluster_sequence : cluster){
        int difference = 0;
        for(int i = 0; i < (int)sequence.size(); ++i){
            if(sequence[i] != cluster_sequence[i]){
                difference++;
                if (difference >= max_cluster_difference){
                    return false;    
                }
            }
        }
    }
    return true;
}

double GVFinder::calculate_average_distance(const std::vector<std::string> &cluster, const std::string &sequence, const int &max_cluster_difference){
    double sum = 0;
    for(const auto &cluster_sequence : cluster){
        int difference = 0;
        for(int i = 0; i < (int)sequence.size(); ++i){
            if(sequence[i] != cluster_sequence[i]){
                difference ++;
                if (difference >= max_cluster_difference){
                    return -1;    
                }
            }
        }
        sum += difference;
    }
    return sum / cluster.size();
}




void GVFinder::cluster_msa_1(const std::vector<std::string> &msa){
    for(const auto &sequence : msa) {
        bool found_cluster = false;
        for(auto &cluster : clusters) {
            if(belongs_to_cluster(cluster, sequence, max_cluster_difference)){
                cluster.push_back(sequence);
                found_cluster = true;
                break;
            }
        }
        if(!found_cluster){
            clusters.push_back({sequence});
        }
    }
}

void GVFinder::cluster_msa_2(const std::vector<std::string> &msa){
    int sequence_index = 0;
    for(const auto &sequence : msa) {
        bool found_cluster = false;
        int cluster_index = 0;
        for(auto &cluster : clusters) {
            if(belongs_to_cluster(cluster, sequence, max_cluster_difference)){
                cluster.push_back(sequence);
                clusters2[cluster_index].push_back(sequences[sequence_index]);
                found_cluster = true;
                break;
            }
            cluster_index++;
        }
        if(!found_cluster){
            clusters.push_back({sequence});
            clusters2.push_back({sequences[sequence_index]});
        }
        sequence_index++;
    }
    clusters = clusters2;
}

void GVFinder::cluster_msa_3(const std::vector<std::string> &msa){
    int sequence_index = 0;
    for(const auto &sequence : msa) {
        bool found_cluster = false;
        double min_avg_distance = max_cluster_difference;
        double current_avg_distance;
        int min_cluster_index = 0;
        
        int cluster_index = 0;
        for(auto &cluster : clusters) {
            current_avg_distance = calculate_average_distance(cluster, sequence, max_cluster_difference);
            if (current_avg_distance != -1){
                if (current_avg_distance < min_avg_distance){
                    min_avg_distance = current_avg_distance;
                    min_cluster_index = cluster_index;
                    found_cluster = true;
                }
            }
            cluster_index++;
        }
        if(!found_cluster){
            clusters.push_back({sequence});
            clusters2.push_back({sequences[sequence_index]});
        } else {
            clusters[min_cluster_index].push_back(sequence);
            clusters2[min_cluster_index].push_back(sequences[sequence_index]);
        }
        sequence_index++;
    }
    clusters = clusters2;
}

std::string GVFinder::get_consensus(const std::vector<std::string> &cluster){
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 0, -1, method == 1 ? -100 : -1);  // linear gaps match mismatch gap
    spoa::Graph graph{};

    for (const auto& it : cluster) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    return graph.GenerateConsensus();
}

void GVFinder::calculate_results(){
    std::sort(clusters.begin(), clusters.end(), compare_by_size);
    while(clusters.back().size() < 10)
        clusters.pop_back();

    for(const auto &cluster : clusters){
        results.push_back(get_consensus(cluster));
    }
}

int GVFinder::get_max_difference(const std::string &X, const std::string &Y){
    if(X.size() != Y.size())
        return -1;
    
    int difference = 0;
    for(int i = 0; i < (int)X.size(); ++i)
        if(X[i] != Y[i])
            difference ++;
    return difference;
}

std::vector<std::string> GVFinder::get_known_results(std::string path) {
    std::vector<std::string> results;
    std::ifstream stream;

    stream.open(path);
    std::string line;
    while(std::getline(stream, line)){
        //if(std::getline(stream, line)){
            results.push_back(line);
        //}
    }
    stream.close();
    return results;
}

void GVFinder::compare_with_known_results(const std::vector<std::string> &known_results, 
                                        const std::vector<std::string> &results) {

    for(const auto &known_result : known_results){
        std::cout << "Known: " << known_result << "\n";
    }

    for(const auto &known_result : known_results){
        int min_difference = 1e9;
        for(const auto &result : results){
            min_difference = std::min(min_difference, get_max_difference(known_result, result));
        }
        std::cout << "Min difference: " << min_difference << "\n";
    }
}

void GVFinder::output_to_file(const std::vector<std::string> &results, std::string path) {
    std::ofstream stream;

    stream.open(path);
    for(const auto &result : results) {
        stream << result << "\n";
    }
    stream.close();
}

void GVFinder::output(std::string path="") {
    std::cerr << "Clusters size: " << clusters.size() << "\n";
    for(int i = 0; i < (int)results.size(); ++i){
        std::cout << "Size: " << clusters[i].size() << "\n";
        std::cout << "Consensus: " << results[i] << "\n";
    }

    if(path != "")
        output_to_file(results, path);

    //std::vector<std::string> known_results = get_known_results(path);
    //compare_with_known_results(known_results, results);
    //std::cout << "Usporedba izmedu sebe:\n";
    //compare_with_known_results(results, results);
}

void GVFinder::solve() {
    find_alignment();
    calculate_results();
}

/*void GVFinder::compare_results(std::string path1, std::string path2){

}*/