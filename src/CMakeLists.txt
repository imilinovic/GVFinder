add_executable(${PROJECT_NAME} 
    main.cpp
    GVFinder.cpp
    ${CMAKE_SOURCE_DIR}/tests/test_helper.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
)

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wpedantic -Werror)

target_link_libraries(${PROJECT_NAME} 
    spoa::spoa
)