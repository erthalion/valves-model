cmake_minimum_required(VERSION 2.6)
set(BIN ValveModel)
project(${BIN} CXX)
set(SOURCE_ROOT src)
set(BUILD_ROOT build)
aux_source_directory(${SOURCE_ROOT} SOURCES)

link_directories("lib")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-o0 -g -DDEBUG -Wno-unused-result")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-o3 -march=native -mtune=native -funsafe-math-optimizations -fforce-addr -fomit-frame-pointer -foptimize-sibling-calls -fno-stack-protector -finline-small-functions -fprefetch-loop-arrays -funsafe-loop-optimizations -funroll-loops -Wno-unused-result")
endif(CMAKE_BUILD_TYPE STREQUAL "Release")

if(CMAKE_PROFILE_TYPE STREQUAL "Generate")
    set(CMAKE_EXE_LINKER_FLAGS "-fprofile-generate")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
endif(CMAKE_PROFILE_TYPE STREQUAL "Generate")

if(CMAKE_PROFILE_TYPE STREQUAL "Use")
    set(CMAKE_EXE_LINKER_FLAGS "-fprofile-use")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-use")
endif(CMAKE_PROFILE_TYPE STREQUAL "Use")

add_executable( ${BUILD_ROOT}/${BIN}.${CMAKE_BUILD_TYPE}  ${SOURCES} )
target_link_libraries (${BUILD_ROOT}/${BIN}.${CMAKE_BUILD_TYPE} inih)
