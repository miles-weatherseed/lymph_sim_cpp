cmake_minimum_required(VERSION 3.15)
project(lymph_sim)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS -pthread)

add_executable(lymph_sim simulation.cpp classesAndFunctions.h simulation.h parametersAndLibraries.h csv.h)