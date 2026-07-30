cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents_davies_85_parallel_variable_parameters LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents_davies_85_parallel_variable_parameters main.cpp 
    solver.cpp
)

add_custom_target(copy_opencl ALL
     COMMAND ${CMAKE_COMMAND} -E copy_directory
             ${CMAKE_SOURCE_DIR}/opencl
             ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/opencl
     DEPENDS ${CMAKE_SOURCE_DIR}/opencl
     COMMENT "Copying opencl into binary directory"
 )

add_dependencies(wind_induced_currents_davies_85_parallel_variable_parameters copy_opencl)

target_link_libraries(wind_induced_currents_davies_85_parallel_variable_parameters PRIVATE
    OpenCL::OpenCL
)
