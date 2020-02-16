# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("derivative")

include_directories(derivative/include)

file(GLOB SOURCES "derivative/*.cpp")
add_executable(test ${SOURCES})

