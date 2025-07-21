cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents_davies_85_constant_parameters LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents_davies_85_constant_parameters main.cpp 
    solver.cpp
)

target_link_libraries(wind_induced_currents_davies_85_constant_parameters PRIVATE
    slae
)
