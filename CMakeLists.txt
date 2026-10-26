cmake_minimum_required(VERSION 3.11.0)
set(CMAKE_C_CSTANDARD C11)
project(SAMPLE C)
set(CMAKE_C_FLAGS "-O2 -Wall -Werror")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_executable(makefactortable makefactortable.c)
add_executable(sample sample.c)

if(WIN32)
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/win32/mman.o
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/win32
        COMMAND gcc -O2 -o mman.o -c mman.c
    )
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/factor.bin
        DEPENDS ${PROJECT_SOURCE_DIR}/makefactortable
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND curl https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
        COMMAND 7z x factor64_data.xz
        COMMAND makefactortable >> factor64_data
        COMMAND icacls factor_data /grant *S-1-1-0:R
        COMMAND move factor64_data factor.bin
    )
else()
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/factor.bin
        DEPENDS ${PROJECT_SOURCE_DIR}/makefactortable
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND wget https://people.maths.bris.ac.uk/~maarb/code/factor64_data.xz
        COMMAND unxz factor64_data.xz
        COMMAND makefactortable >> factor64_data
        COMMAND chmod 444 factor64_data
        COMMAND mv factor64_data factor.bin
    )
endif()

add_custom_target(
    FACTOR_BIN
    DEPENDS ${PROJECT_SOURCE_DIR}/factor.bin
)

add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/factor64.o
    DEPENDS ${PROJECT_SOURCE_DIR}/factor64.c
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND gcc -O2 -march=native -c ${PROJECT_SOURCE_DIR}/factor64.c
)

add_custom_target(
    FACTOR64_O
    DEPENDS ${PROJECT_SOURCE_DIR}/factor64.o FACTOR_BIN
)

add_library(libfactor64 OBJECT IMPORTED GLOBAL)

if (WIN32)
    add_custom_target(
        MMAN_O
        DEPENDS ${PROJECT_SOURCE_DIR}/win32/mman.o
    )
    add_dependencies(FACTOR64_O MMAN_O)
    set_property(TARGET libfactor64 PROPERTY IMPORTED_OBJECTS ${PROJECT_SOURCE_DIR}/factor64.o ${PROJECT_SOURCE_DIR}/win32/mman.o)
else()
    set_property(TARGET libfactor64 PROPERTY IMPORTED_OBJECTS ${PROJECT_SOURCE_DIR}/factor64.o)
endif()

add_dependencies(sample FACTOR64_O)

target_link_libraries(sample libfactor64)
