cmake_minimum_required(VERSION 3.16)

project(laplace LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(laplace main.cpp)

target_link_libraries(laplace PRIVATE
    slae
)
