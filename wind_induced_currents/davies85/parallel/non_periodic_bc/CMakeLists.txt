cmake_minimum_required(VERSION 3.16)

project(wind_induced_currents_davies_85_parallel_non_periodic_bc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wind_induced_currents_davies_85_parallel_non_periodic_bc 
    main.cpp
    solver.cpp

    writers/height.cpp
    writers/surface.cpp
    writers/volume.cpp
    writers/statistics.cpp

    generators/area/uniform_generator.cpp
    generators/area/proj_generator.cpp

    generators/dz/uniform_generator.cpp
    generators/dz/triple_point_generator.cpp

    generators/bathymetry/uniform_generator.cpp
    generators/bathymetry/gebco_generator.cpp

    generators/wind/uniform_generator.cpp
    generators/wind/ecmwf_generator.cpp

    generators/viscosity/uniform_generator.cpp
    generators/viscosity/wind_speed_generator.cpp
    generators/viscosity/linear_generator.cpp
)

target_copy_opencl(wind_induced_currents_davies_85_parallel_non_periodic_bc)
target_copy_bathymetry_data(wind_induced_currents_davies_85_parallel_non_periodic_bc)
target_copy_wind_data(wind_induced_currents_davies_85_parallel_non_periodic_bc)

target_finalize_asset_copy(wind_induced_currents_davies_85_parallel_non_periodic_bc)

target_add_manifest(wind_induced_currents_davies_85_parallel_non_periodic_bc)

target_link_libraries(wind_induced_currents_davies_85_parallel_non_periodic_bc PRIVATE
    utils
    OpenCL::OpenCL
    netCDF::netcdf netCDF::netcdf-cxx4
    PROJ::proj
    Boost::math        
)
