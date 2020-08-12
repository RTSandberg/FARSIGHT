cmake_minimum_required(VERSION 3.10)

# set the project name
project(FARRSIGHT)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release" "None")
endif()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(OpenMP REQUIRED)

add_subdirectory(src)
include_directories(inc)
include_directories(inc_eigen/eigen-3.3.7)
link_directories(src)

# add the executable

add_executable(farrsight run_FARRSIGHT.cpp)
target_link_libraries(farrsight AMRObj AMRInterpObj AMRStepObj AMRFieldObj AMRIOObj AMRMeshObj PanelObj)

target_link_libraries(farrsight OpenMP::OpenMP_CXX)




install(TARGETS farrsight DESTINATION bin)


enable_testing()

add_test(NAME Runs COMMAND farrsight 1 0.0 1.0 -2.0 2.0)