add_executable(${PROJECT_NAME} 
    main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
)

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -Wextra -Wpedantic -Werror)

target_link_libraries(${PROJECT_NAME} 
    spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>
    spoa::spoa
)