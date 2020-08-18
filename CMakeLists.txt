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

option(ENABLE_GPU_BUILD "build OpenACC GPU version"             OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(OpenMP REQUIRED)

# from Leighton
#Some PGI versions seem to require a link flag for OpenMP
set_property(TARGET OpenMP::OpenMP_C APPEND
                                     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         $<$<CXX_COMPILER_ID:PGI>:-mp>) 
#Creating an interface for OpenACC
add_library(OpenACC_CXX INTERFACE)
set_property(TARGET OpenACC_CXX     PROPERTY
                                         INTERFACE_COMPILE_OPTIONS
                                         -acc
#                                         -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
                                         -Minfo)
set_property(TARGET OpenACC_CXX     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         -acc
#                                         -ta=tesla:ccall,cc60
                                        )

if(ENABLE_GPU_BUILD AND NOT (CMAKE_C_COMPILER_ID MATCHES PGI))
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "build OpenACC GPU version" FORCE)
endif()


add_subdirectory(src)
include_directories(inc)
include_directories(inc_eigen/eigen-3.3.7)
link_directories(src)

# add the executable
set(EXE farrsight)
add_executable(${EXE} run_FARRSIGHT.cpp)
target_link_libraries(${EXE}    PanelObj 
                                icsObj)
if (ENABLE_GPU_BUILD)
    target_link_libraries(farrsight AMRObj_gpu)
else()
    target_link_libraries(farrsight AMRObj_cpu)
endif()

target_link_libraries(farrsight OpenMP::OpenMP_CXX)

install(TARGETS farrsight DESTINATION bin)


enable_testing()

add_test(NAME Runs COMMAND farrsight 1 0.0 1.0 -2.0 2.0)