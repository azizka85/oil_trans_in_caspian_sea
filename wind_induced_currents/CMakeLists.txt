cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents main.cpp 
    davies85/solver.cpp
    davies85/variable_step_size/solver.cpp
)

target_link_libraries(wind_induced_currents PRIVATE
    slae
)
