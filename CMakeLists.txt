cmake_minimum_required(VERSION 3.16)
project(Hedge_cpp)

set(CMAKE_CXX_STANDARD 14)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

add_executable(Hedge_cpp main.cpp src/utils.h src/utils.cpp)
