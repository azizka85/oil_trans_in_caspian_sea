cmake_minimum_required(VERSION 3.16)

project(oil_trans_in_caspian_sea LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(slae)

add_subdirectory(wind_induced_currents/davies85)
