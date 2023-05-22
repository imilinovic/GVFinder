#include <gtest/gtest.h>
#include "test_helper.hpp"

Test_helper helper;

TEST(test_big_punishment, test_sample_29)
{
    Test_result result = helper.solve("data/J29_B_CE_IonXpress_005.fastq", 2, 29);
    std::cout << "Number of clusters: " << result.number_of_clusters << "\n";
    ASSERT_EQ(result.number_of_matches, 2);
}

TEST(test_big_punishment, test_sample_30)
{
    Test_result result = helper.solve("data/J30_B_CE_IonXpress_006.fastq", 2, 30);
    std::cout << "Number of clusters: " << result.number_of_clusters << "\n";
    ASSERT_EQ(result.number_of_matches, 2);
}

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}