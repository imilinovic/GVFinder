#pragma once

struct Test_result{
    int number_of_matches;
    int number_of_clusters;
    Test_result(){};
    Test_result(int number_of_matches, int number_of_clusters) : number_of_matches(number_of_matches), number_of_clusters(number_of_clusters) {};
};