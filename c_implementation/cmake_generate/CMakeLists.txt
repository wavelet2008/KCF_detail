cmake_minimum_required(VERSION 2.8)
project(kcf_without_opencv)

if(NOT WIN32)
ADD_DEFINITIONS("-std=c++0x -O3")
endif(NOT WIN32)

include_directories(src)
FILE(GLOB_RECURSE sourcefiles "../src/*.c" "../src/*.h")
add_executable( KCF ${sourcefiles} )
target_link_libraries( KCF )
