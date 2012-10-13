cmake_minimum_required(VERSION 2.6)
set(BIN ValveModel)
project(${BIN} CXX)
set(SOURCE_ROOT src)
set(BUILD_ROOT build)
aux_source_directory(${SOURCE_ROOT} SOURCES)

if(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_CXX_FLAGS "-o0 -g")
endif(CMAKE_BUILD_TYPE Debug)

if(CMAKE_BUILD_TYPE Release)
    set(CMAKE_CXX_FLAGS "-o3 -march='athlon' -funsafe-math-optimizations -fforce-addr -fomit-frame-pointer -foptimize-sibling-calls -fno-stack-protector -finline-small-functions -ftree-parallelize-loops=2 -fprefetch-loop-arrays -funsafe-loop-optimizations -floop-parallelize-all -ftree-loop-linear -funroll-loops -ftree-loop-im -fopenmp -std=c++0x -mpreferred-stack-boundary=2")
endif(CMAKE_BUILD_TYPE Release)

add_executable( ${BUILD_ROOT}/${BIN}.${CMAKE_BUILD_TYPE}  ${SOURCES} )
