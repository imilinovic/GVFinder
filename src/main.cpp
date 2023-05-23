#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include "GVFinder.hpp"
#include "spoa/spoa.hpp"
#include "../tests/test_helper.hpp"

int main() {
  // Test_helper helper;
  // Test_result result1, result2;
  // std::ofstream stream;
  // std::string jelenref01 = "GAGCATCTTAAGGCCGAGTGTCATTTCTTCAACGGGACGGAGCGGATGCAGTTCCTGGCGAGATACTTCTATAACGGAGAAGAGTACGCGCGCTTCGACAGCGACGTGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGGACGCCAAGTACTGGAACAGCCAGAAGGAGATCCTGGAGCAGCACCGGGCAGAGGTGGACAGGTACTGCAGACACAACTACGGGGTCGGTGAGAGTTTCACTGTG";
  
  // for(int max_cluster_difference = 5; max_cluster_difference <= 25; max_cluster_difference++){
  //   if (helper.found_sequence("data/J30_B_CE_IonXpress_006.fastq", jelenref01, 1, max_cluster_difference, 0)){
  //     std::cout << "Uspjeh\n";
  //     stream.open("combo.txt", std::ios_base::app);
  //     stream << 1 << ' ' << max_cluster_difference << ' ' << 0 << ' ' << "\n";
  //     stream.close();
  //   }
  //   for(int max_size_difference = 0; max_size_difference <= 15; max_size_difference++){
  //     for(int method = 2; method <= 3; method++){
  //       if (helper.found_sequence("data/J30_B_CE_IonXpress_006.fastq", jelenref01, method, max_cluster_difference, max_size_difference)){
  //         std::cout << "Uspjeh\n";
  //         stream.open("combo.txt", std::ios_base::app);
  //         stream << method << ' ' << max_cluster_difference << ' ' << max_size_difference << ' ' << "\n";
  //         stream.close();
  //       }
  //      }
  //     }
  //   }
  

  // Test_helper helper;
  // Test_result result1, result2;
  // std::ofstream stream;

  // for(int max_cluster_difference = 5; max_cluster_difference <= 25; max_cluster_difference++){
  //   result1 = helper.solve("data/J29_B_CE_IonXpress_005.fastq", 1, 29, false, max_cluster_difference, 0);
  //   result2 = helper.solve("data/J30_B_CE_IonXpress_006.fastq", 1, 30, false, max_cluster_difference, 0);
  //   if(result1.number_of_matches >= 2 && result2.number_of_matches >= 2){
  //     std::cout << "Uspjeh\n";
  //     stream.open("combo.txt", std::ios_base::app);
  //     stream << 1 << ' ' << max_cluster_difference << ' ' << 0 << ' '
  //       << result1.number_of_clusters << ' ' << result2.number_of_clusters << ' '
  //       << result1.number_of_matches << ' ' << result2.number_of_matches << "\n";
  //     stream.close();
  //   }
  //   for(int max_size_difference = 0; max_size_difference <= 10; max_size_difference++){
  //     for(int method = 2; method <= 3; method++){
  //       result1 = helper.solve("data/J29_B_CE_IonXpress_005.fastq", method, 29, false, max_cluster_difference, max_size_difference);
  //       result2 = helper.solve("data/J30_B_CE_IonXpress_006.fastq", method, 30, false, max_cluster_difference, max_size_difference);
  //       if(result1.number_of_matches >= 2 && result2.number_of_matches >= 2){
  //         std::cout << "Uspjeh\n";
  //         stream.open("combo.txt", std::ios_base::app);
  //         stream << method << ' ' << max_cluster_difference << ' ' << max_size_difference << ' '
  //           << result1.number_of_clusters << ' ' << result2.number_of_clusters << ' '
  //           << result1.number_of_matches << ' ' << result2.number_of_matches << "\n";
  //         stream.close();
  //       }
  //     }
  //   }
  // }

  std::string line;
  
  int method;
  int max_size_difference;
  int max_cluster_difference;
  std::cin >> method;
  std::cin >> max_size_difference;
  std::cin >> max_cluster_difference;
    GVFinder gvfinder = GVFinder("data/J30_B_CE_IonXpress_006.fastq", method, max_size_difference, max_cluster_difference);
    // std::cout << gvfinder.get_allele_size() << " " << gvfinder.get_data_size() << "\n";
    gvfinder.solve();
    gvfinder.output("");
    gvfinder.check("GAGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCGGGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGATTTCATGGAGCAGAAGCGGGCCGAGGTGGACACGGTGTGCAGACACAACTACGGGGTTATTGAGAGTTTCACTGTG", "jelenref02");
    gvfinder.check("GAGCATCTTAAGGCCGAGTGTCATTTCTTCAACGGGACGGAGCGGATGCAGTTCCTGGCGAGATACTTCTATAACGGAGAAGAGTACGCGCGCTTCGACAGCGACGTGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGGACGCCAAGTACTGGAACAGCCAGAAGGAGATCCTGGAGCAGCACCGGGCAGAGGTGGACAGGTACTGCAGACACAACTACGGGGTCGGTGAGAGTTTCACTGTG", "jelenref01");
    gvfinder.check("ATGTATACTAAGAAAGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGCTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTTCCGGGCGGTGACCGAGCTGGGGCGGCCGGCGGCCGAGGGCTGGAACAGCCAGAAGGAGCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACGTACTGCAGACACAACTACGGGGTTATTGAGAGTTTCACTGTG", "jelenref04");

    gvfinder = GVFinder("data/J29_B_CE_IonXpress_005.fastq", method, max_size_difference, max_cluster_difference);
    // std::cout << gvfinder.get_allele_size() << " " << gvfinder.get_data_size() << "\n";
    gvfinder.solve();
    gvfinder.output("");
    gvfinder.check("CTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTG", "jelenref05");
    gvfinder.check("CTGTATACTACGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTGGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACAGGTACTGCAGACACAACTACGGGGTTCTTGACAGTTTCGCTGTG", "jelenref06");
    gvfinder.check("GAGCATCATAAGTGCGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCAGTTCCTGCAGAGATACATCTATAACCGGGAAGAGTACGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGTCCGCCAAGTACTATAACAGCCAGAAGGAGCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAGGTACTGCAGACACAACTACGGGGTCGTTGAGAGTTTCACTGTG", "jelenref07");




  return 0;
}