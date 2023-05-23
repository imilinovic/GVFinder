enable_testing()

add_executable( test_big_punishment
  test_big_punishment.cpp
  ${CMAKE_SOURCE_DIR}/src/GVFinder.cpp
  test_helper.cpp
  utilities.hpp
)

target_link_libraries( test_big_punishment
  GTest::gtest_main
  spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
  spoa::spoa
)

add_executable( test_rough_method
  test_rough_method.cpp
  ${CMAKE_SOURCE_DIR}/src/GVFinder.cpp
  test_helper.cpp
  utilities.hpp
)

target_link_libraries( test_rough_method
  GTest::gtest_main
  spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
  spoa::spoa
)

add_executable( test_closest_clusters
  test_closest_clusters.cpp
  ${CMAKE_SOURCE_DIR}/src/GVFinder.cpp
  test_helper.cpp
  utilities.hpp
)

target_link_libraries( test_closest_clusters
  GTest::gtest_main
  spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
  spoa::spoa
)

include(GoogleTest)
gtest_discover_tests(test_big_punishment)
gtest_discover_tests(test_rough_method)
gtest_discover_tests(test_closest_clusters)