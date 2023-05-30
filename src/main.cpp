#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>
#include "GVFinder.hpp"
#include "spoa/spoa.hpp"
#include "../tests/test_helper.hpp"

int main() {
    std::string filename, output_file;
    std::cout << "Unesite datoteku s podacima: ";
    std::cin >> filename;
    std::cout << "Unesite datoteku za ispis rezultata: ";
    std::cin >> output_file;

    GVFinder gvfinder(filename, 1, 3, 0); 
    gvfinder.solve();
    gvfinder.output(filename, output_file);
}