cmake_minimum_required(VERSION 3.10)
project(TopoMinCut)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(phat)

set(TOPOMINCUT_EIGEN_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/dependencies")

find_package(Boost REQUIRED COMPONENTS graph)

set(TOPOMINCUT_SOURCES
    src/boundary_matrix.cpp
    src/tags.cpp
    src/persistence.cpp
    src/thinning.cpp
    src/iset.cpp
    src/graph_flow.cpp
    src/topomincut_runner.cpp
)

add_library(topomincut_core STATIC ${TOPOMINCUT_SOURCES})

target_include_directories(topomincut_core
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${TOPOMINCUT_EIGEN_INCLUDE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/phat/include"
        ${Boost_INCLUDE_DIRS}
)

target_link_libraries(topomincut_core PUBLIC phat ${Boost_LIBRARIES})

add_executable(TopoMinCut src/topo_min_cut_cli.cpp)
target_link_libraries(TopoMinCut PRIVATE topomincut_core)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(topomincut_core PRIVATE -Wall -Wextra -O2)
    target_compile_options(TopoMinCut PRIVATE -Wall -Wextra -O2)
endif()
