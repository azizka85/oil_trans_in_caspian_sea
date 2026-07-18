cmake_minimum_required(VERSION 3.16)

project(diffusion_eulerian_2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(diffusion_eulerian_2d main.cpp 
    solver.cpp
)
