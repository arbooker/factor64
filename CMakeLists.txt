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
