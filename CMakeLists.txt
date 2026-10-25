cmake_minimum_required(VERSION 3.11.0)
set(CMAKE_C_CSTANDARD C11)
project(SAMPLE C)
set(CMAKE_C_FLAGS "-O2 -Wall -Werror")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(makefactortable makefactortable.c)

if(WIN32)
else()
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/factor.bin
        DEPENDS ${PROJECT_SOURCE_DIR}/makefactortable
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND wget https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
        COMMAND unxz factor64_data.xz
        COMMAND ${PROJECT_SOURCE_DIR}/bin/makefactortable >> factor64_data
        COMMAND chmod 444 factor64_data
        COMMAND mv factor64_data factor.bin
    )
endif()

add_custom_target(
    FACTOR_BIN
    DEPENDS ${PROJECT_SOURCE_DIR}/factor.bin
)

if(WIN32)
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/factor64.o
        DEPENDS ${PROJECT_SOURCE_DIR}/factor64.c ${PROJECT_SOURCE_DIR}/win32/mman.c ${PROJECT_SOURCE_DIR}/win32/mman.h
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND gcc -O2 -march=native -c ${PROJECT_SOURCE_DIR}/factor64.c ${PROJECT_SOURCE_DIR}/win32/mman.c
    )
else()
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/factor64.o
        DEPENDS ${PROJECT_SOURCE_DIR}/factor64.c
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND gcc -O2 -march=native -c ${PROJECT_SOURCE_DIR}/factor64.c
    )
endif()

add_custom_target(
    FACTOR64_O
    DEPENDS ${PROJECT_SOURCE_DIR}/factor64.o FACTOR_BIN
)

add_executable(sample sample.c)
add_dependencies(sample FACTOR64_O)

add_library(libfactor64 OBJECT IMPORTED GLOBAL)
set_property(TARGET libfactor64 PROPERTY IMPORTED_OBJECTS ${PROJECT_SOURCE_DIR}/factor64.o)
target_link_libraries(sample libfactor64)
