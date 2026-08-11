cmake_minimum_required(VERSION 3.16)

project(oil_trans_in_caspian_sea LANGUAGES CXX)

find_package(OpenCL CONFIG REQUIRED)
find_package(netCDFCxx CONFIG REQUIRED)
find_package(PROJ CONFIG REQUIRED)
find_package(boost_math CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_file(PROJ_DB_FILE_PATH 
    NAMES proj.db
    PATHS "${CMAKE_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET}"
    PATH_SUFFIXES "share/proj"
    NO_DEFAULT_PATH
    REQUIRED
)

function(target_copy_opencl TARGET_NAME)
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY ASSETS_TO_COPY "opencl")
endfunction()

function(target_copy_bathymetry_data TARGET_NAME)
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY ASSETS_TO_COPY "data/bathymetry")
endfunction()

function(target_copy_wind_data TARGET_NAME)
    set_property(TARGET ${TARGET_NAME} APPEND PROPERTY ASSETS_TO_COPY "data/wind")
endfunction()

function(target_finalize_asset_copy TARGET_NAME)
    get_target_property(ASSET_LIST ${TARGET_NAME} ASSETS_TO_COPY)
    
    if(ASSET_LIST)
        set(COPY_COMMANDS 
            COMMAND ${CMAKE_COMMAND} -E copy
                    "${PROJ_DB_FILE_PATH}"
                    "$<TARGET_FILE_DIR:${TARGET_NAME}>/proj.db"
        )    

        foreach(ASSET_DIR ${ASSET_LIST})
            list(APPEND COPY_COMMANDS 
                 COMMAND ${CMAKE_COMMAND} -E copy_directory
                         "${CMAKE_SOURCE_DIR}/${ASSET_DIR}"
                         "$<TARGET_FILE_DIR:${TARGET_NAME}>/${ASSET_DIR}"
            )
        endforeach()

        add_custom_command(
            TARGET ${TARGET_NAME}
            POST_BUILD
            ${COPY_COMMANDS}
            COMMENT "Copying registered assets for ${TARGET_NAME}"
        )
    endif()
endfunction()

function(target_add_manifest TARGET_NAME)
    if(MSVC)
        target_sources(${TARGET_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/app.manifest")
    endif()
endfunction()

add_subdirectory(utils)

add_subdirectory(slae)

add_subdirectory(wind_induced_currents/davies85/constant_parameters)
add_subdirectory(wind_induced_currents/davies85/variable_parameters)
add_subdirectory(wind_induced_currents/davies85/non_periodic_bc)
add_subdirectory(wind_induced_currents/davies85/parallel/variable_parameters)
add_subdirectory(wind_induced_currents/davies85/parallel/non_periodic_bc)

add_subdirectory(diffusion/eulerian/1d)
add_subdirectory(diffusion/eulerian/2d)

add_subdirectory(diffusion/lagrangian/1d)
