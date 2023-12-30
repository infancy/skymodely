cmake_minimum_required(VERSION 3.15)

project(sky)

set(CMAKE_CXX_STANDARD 23)

add_executable(skymodely skymodely.cpp Arhosek12/ArHosekSkyModel.cpp)
target_compile_options(skymodely PUBLIC -fopenmp)
