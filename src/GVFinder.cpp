#include "GVFinder.hpp"
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include "spoa/spoa.hpp"

GVFinder::GVFinder(std::string data_path) {
    std::vector<std::string> data;
    std::ifstream stream;
    std::unordered_map<int, int> length_cnt;

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
        //std::cout << len.first << ": " << len.second << "\n";
        if(len.second > max_len_cnt){
            max_len_cnt = len.second;
            most_common_len = len.first;
        }
    }

    for(auto allele : data){
        if((int)allele.size() == most_common_len)
            sequences.push_back(allele);
    }
}

void GVFinder::find_alignment() {
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 0, -1, -100);  // linear gaps match mismatch gap
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
    
    cluster_msa(msa);
    /*for (const auto& it : msa) {
        std::cerr << it << "  VELICINA: " << it.size() << std::endl;
    }*/
}

bool GVFinder::compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y){
    return X.size() > Y.size();
}

int GVFinder::get_max_difference(const std::vector<std::string> &cluster, const std::string &sequence){
    int max_difference = 0;
    for(const auto &cluster_sequence : cluster){
        int difference = 0;
        for(int i = 0; i < (int)sequence.size(); ++i)
            if(sequence[i] != cluster_sequence[i])
                difference ++;
        max_difference = std::max(max_difference, difference);
    }
    return max_difference;
}

void GVFinder::cluster_msa(const std::vector<std::string> &msa){
    for(const auto &sequence : msa) {
        bool found_cluster = false;
        for(auto &cluster : clusters) {
            int max_difference = get_max_difference(cluster, sequence);
            if(max_difference < max_cluster_difference){
                cluster.push_back(sequence);
                found_cluster = true;
                break;
            }
        }
        if(!found_cluster){
            clusters.push_back({sequence});
        }
    }
    std::sort(clusters.begin(), clusters.end(), compare_by_size);

    std::cerr << "Clusters size: " << clusters.size() << "\n";
    for(const auto &cluster : clusters)
        std::cout << "Size: " << cluster.size() << "\n";
}