enable_testing()

add_executable( test_1
  test_1.cpp
)

target_link_libraries( test_1
  GTest::gtest_main
  spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
)

include(GoogleTest)
gtest_discover_tests(test_1)