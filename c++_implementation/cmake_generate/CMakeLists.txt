cmake_minimum_required(VERSION 2.8)
project(test)

set(OpenCV_DIR  C:/opencv/build/x64/vc14/lib)
find_package(OpenCV REQUIRED)

if(NOT WIN32)
ADD_DEFINITIONS("-std=c++0x -O3")
endif(NOT WIN32)

include_directories(src)
FILE(GLOB_RECURSE sourcefiles "../src/*.cpp" "../src/*.hpp" "../src/*.h" )
add_executable( KCF ${sourcefiles} )
target_link_libraries( KCF ${OpenCV_LIBS})
