cmake_minimum_required(VERSION 3.25)
project(untitled53 C)

set(CMAKE_C_STANDARD 11)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

add_executable(untitled53 main.c matrix.h matrix.c)
