cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents_davies_85 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents_davies_85 main.cpp 
    solver.cpp
    variable_step_size/solver.cpp
    variable_step_size/variable_depth/solver.cpp
)

target_link_libraries(wind_induced_currents_davies_85 PRIVATE
    slae
)
