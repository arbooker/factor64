cmake_minimum_required(VERSION 3.11.0)
set(CMAKE_C_CSTANDARD C11)
project(SAMPLE C)
set(CMAKE_C_FLAGS "-O2 -Wall -Werror")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/lib/factor64.bin
    DEPENDS ${PROJECT_SOURCE_DIR}/bin/makefactortable
    COMMAND wget -nc -P ${PROJECT_SOURCE_DIR}/factor64_data https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
    COMMAND unxz -k ${PROJECT_SOURCE_DIR}/factor64_data/factor64_data.xz
    COMMAND ${PROJECT_SOURCE_DIR}/bin/makefactortable >> ${PROJECT_SOURCE_DIR}/factor64_data/factor64_data
    COMMAND chmod 444 ${PROJECT_SOURCE_DIR}/factor64_data/factor64_data
    COMMAND mv ${PROJECT_SOURCE_DIR}/factor64_data/factor64_data ${PROJECT_SOURCE_DIR}/lib/factor.bin
)

add_executable(makefactortable makefactortable.c)
add_executable(sample sample.c)
add_library(factor64 SHARED factor64.c ${PROJECT_SOURCE_DIR}/lib/factor64.bin)
target_link_libraries(sample PRIVATE factor64)
