#[[
    GCC and Cmake are required to compile sample.c
    7zip is also required on Windows. Remember to add it to the path!
    To compile on Linux:
        1. Make a `build` directory in the project root and open a terminal there
        2. Run `cmake ..`
        3. Run `cmake --build .`
    To compile on Windows:
        0. Install GCC using this guide: https://code.visualstudio.com/docs/cpp/config-mingw
        1. Make a `build` directory in the project root and open a terminal there
        2. Run `cmake .. --preset windows`
        3. Run `cmake --build .`
]]

cmake_minimum_required(VERSION 3.11.0)
set(CMAKE_C_CSTANDARD C11)
project(SAMPLE C)
set(CMAKE_C_FLAGS "-O2 -Wall -Werror")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(makefactortable makefactortable.c)
add_executable(sample sample.c)
add_library(factor64 factor64.c)
if(WIN32)
    target_sources(factor64 PRIVATE win32/mman.c)
endif()

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/factor.bin
    DEPENDS makefactortable
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} -P factor.bin.cmake
)

add_custom_target(
    FACTOR_BIN
    DEPENDS ${PROJECT_SOURCE_DIR}/factor.bin
)

add_dependencies(factor64 FACTOR_BIN)
target_link_libraries(sample factor64)
