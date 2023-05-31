#include "GVFinder.hpp"
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <filesystem>
#include <algorithm>
#include "spoa/spoa.hpp"

// Klasa GVFinder glavna je klasa algoritma, u konstruktoru joj predajemo put (iz root-a) do datoteke s 
// očitanjima, metodu iz {1, 2, 3}, maksimalnu toleranciju za razliku u duljini očitanja (primjenjuje se
// samo za metode 2 i 3) te maksimalnu udaljenost sekvenci kod klasteriranja

GVFinder::GVFinder(std::string data_path, int method, int max_size_difference, int max_cluster_difference) {
    std::vector<std::string> data;
    std::ifstream stream;
    std::unordered_map<int, int> length_cnt;
    this->method = method;
    this->max_size_difference = max_size_difference;
    this->max_cluster_difference = max_cluster_difference;
    
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

// spremamo duljine očitanja u datoteci u rječnik
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

// pronalazimo najčešću duljinu očitanja
    int most_common_len = -1, max_len_cnt = 0;
    for(auto len : length_cnt){
        if(len.second > max_len_cnt){
            max_len_cnt = len.second;
            most_common_len = len.first;
        }
    }

// ako koristimo prvu metodu, koristit ćemo samo očitanja najčešće duljine, a
// inače koristimo očitanja duljine [najčešća - max_size_difference, najčešća + max_size_difference]
    if(method == 1){
        for(auto allele : data){
            if((int)allele.size() == most_common_len){
                sequences.push_back(allele);
            }
        }
    } else{
        for(auto allele : data){    
            if((int)allele.size() <= most_common_len + max_size_difference && (int)allele.size() >= most_common_len - max_size_difference){
                sequences.push_back(allele);
            }
        }
    }
}

// ova funkcija pronalazi višestruko poravnanje odabranih sekvenci te generira klastere
void GVFinder::find_alignment() {
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 0, -1, method == 1 ? -100 : -1);  // parametri korišteni za poravnanje su redom match mismatch, gap
    spoa::Graph graph{};
    // ako koristimo metodu 1, onda se umetanje kažnjava dodatno kako bi sva poravnanja bila jednake(najčešće) duljine

    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    auto consensus = graph.GenerateConsensus();
    auto msa = graph.GenerateMultipleSequenceAlignment();

// ovisno o metodu pozivamo funkciju za klasteranje
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

// funkcija uspoređuje vektore stringova po veličinama 
bool GVFinder::compare_by_size(const std::vector<std::string> &X, const std::vector<std::string> &Y){
    return X.size() > Y.size();
}

// funkcija vraća true ako je string sequence od svakog stringa u clusteru udaljen manje od max_cluster_difference
// udaljenost se računa kao broj pozicija s različitim znakovima
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

//funkcija se koristi u trećoj metodi i vraća prosječnu udaljenost sekvence od svih sekvenci u klasteru,
// te vraća -1 ako sekvenca ne pripada klasteru (udaljena je od neke sekvence za max_cluster_difference ili više)
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

// metoda dodaje sekvencu u prvi klaster kojemu sekvenca pripada ili stvara novi klaster ako sekvenca ne pripada niti jednom
// postojećem klasteru
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

// budući da metoda 2 radi gradi višestruko poravnanje od sekvenci različitih duljina, ta će poravnanja
// sadržavati umetanja i brisanja. Kako bi iz nastalih klastera mogli generirati konsenzusnu sekvencu 
// bez umetanja i brisanja, ova metoda gradi dvije skupine klastera - u jednoj su klasteri s očitanjima iz
// višestrukog poravnanja (sadrže umetanja i brisanja) i njih koristimo za klasteriranje, a u drugoj su 
// originalna očitanja iz kojih su nastale odgovarajuće sekvence iz prve skupine. Drugu skupinu koristimo 
// za generiranje konsenzusa klastera. ostatak algoritma isti je kao kod prvog, sekvenca se dodaje u
// prvi klaster kojemu pripada ili se stvara novi klaster
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

// treća metoda isto kao i druga koristi očitanja različitih duljina pa su joj potrebna dva skupa klastera
// ono u čemu se treća metoda razlikuje od druge je to što ne dodaje sekvencu u prvi klaster kojemu pripada, 
// ju dodaje u klaster čije su sekvence u prosjeku najbliže toj sekvenci (ako mu sekvenca uopće pripada)
// za računanje prosječne udaljenosti koristi se funkcija calculate_average_distance
// ako ne pripada niti jednom klasteru stvara se novi
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

// funkcija prima vektor stringova (klaster) i generira višestruko poravnanje kako bi iz njega
// generirao konsenzus koji će biti predstavnik tog klastera
std::string GVFinder::get_consensus(const std::vector<std::string> &cluster){
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 0, -1, method == 1 ? -100 : -1);  // match mismatch gap
    spoa::Graph graph{};

    for (const auto& it : cluster) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    return graph.GenerateConsensus();
}

// funkcija uklanja klastere koji sadrže 10 ili manje sekvenci te od preostalih bira 
// četiri najveća i sprema ih u result
void GVFinder::calculate_results(){
    std::sort(clusters.begin(), clusters.end(), compare_by_size);
    while(clusters.size() > 4 || clusters.back().size() < 10)
        clusters.pop_back();

    for(const auto &cluster : clusters){
        results.push_back(get_consensus(cluster));
    }
}

// funkcija na kraju nije iskorištena
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
        results.push_back(line);
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

// funkcija zapisuje rezultate u datoteku
void GVFinder::output_to_file(std::string filename, const std::vector<std::string> &results, std::string path) {
    std::ofstream stream;

    stream.open(path, std::ios_base::app);
    stream << filename << "\n";
    for(const auto &result : results) {
        stream << result << "\n";
    }
    stream << "\n";
    stream.close();
}

// funkcija ispisuje rezultate na standardni izlaz
void GVFinder::output(std::string filename, std::string path) {
    for(int i = 0; i < (int)results.size(); ++i){
        std::cout << "Cluster size: " << clusters[i].size() << "\n";
        std::cout << "Cluster consensus: " << results[i] << "\n";
        std::cout << "\n";
    }

    if(path != "")
        output_to_file(filename, results, path);
}

void GVFinder::check(std::string sequence, std::string name) {
    for (int i = 0; i < int(results.size()); i++){
        if (results[i].find(sequence) != std::string::npos){
            std::cout << "found " << name << " in " << i + 1 << ". cluster" << "\n";
            return; 
        }
    }
}

void GVFinder::solve() {
    find_alignment();
    calculate_results();
}