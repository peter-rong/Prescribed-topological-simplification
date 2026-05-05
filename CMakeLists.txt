cmake_minimum_required(VERSION 3.10)
project(PrescribedTopologicalSimplification)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(PrescribedTopologicalSimplification
    main.cpp
    cubical_mode.cpp
    tet_mode.cpp
    marching_cubes.cpp
)

target_include_directories(PrescribedTopologicalSimplification PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

# TopoTet optionally uses parallel std::sort via <execution>; disable if your toolchain lacks TBB linkage.
# target_compile_definitions(PrescribedTopologicalSimplification PRIVATE TOPOTET_NO_PAR_SORT)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(PrescribedTopologicalSimplification PRIVATE -Wall -Wextra -O2)
endif()
