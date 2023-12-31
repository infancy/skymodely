cmake_minimum_required(VERSION 3.15)

project(skymodely)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

add_executable(skymodely skymodely.cpp Arhosek12/ArHosekSkyModel.c)
target_compile_options(skymodely PUBLIC -fopenmp)
