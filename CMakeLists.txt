cmake_minimum_required(VERSION 3.11.0)
set(CMAKE_C_CSTANDARD C11)
project(SAMPLE C)
set(CMAKE_C_FLAGS "-O2 -Wall -Werror")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/factor.bin.d/factor.bin
    DEPENDS ${PROJECT_SOURCE_DIR}/bin/makefactortable
    COMMAND wget https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/factor.bin.d
    COMMAND unxz factor64_data.xz
    COMMAND ${PROJECT_SOURCE_DIR}/bin/makefactortable >> factor64_data
    COMMAND chmod 444 factor64_data
    COMMAND mv factor64_data factor.bin
)

add_custom_target(
    FACTOR_BIN ALL
    DEPENDS ${PROJECT_SOURCE_DIR}/factor.bin.d/factor.bin
)

add_executable(makefactortable makefactortable.c)
add_executable(sample sample.c)
add_library(factor64 SHARED factor64.c)
add_dependencies(factor64 FACTOR_BIN)
target_link_libraries(sample PUBLIC factor64)
