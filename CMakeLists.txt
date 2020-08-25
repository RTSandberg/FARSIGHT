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
find_package(MPI    REQUIRED)

set(CMAKE_PREFIX_PATH "/usr/local")
#option(BaryTree_DIR     /usr/local/)
#set(BaryTree_DIR ${CMAKE_INSTALL_PREFIX})
find_library(BaryTree_lib NAMES BaryTree_cpu
                HINTS /usr/local/lib)
#find_package(BaryTree REQUIRED)
#add_library(BaryTree_Interface INTERFACE)
#target_link_libraries(BaryTree_Interface INTERFACE barytree)
#target_include_directories(BaryTree_Interface INTERFACE ${BaryTree_INCLUDE_DIRS})
#    find_package(Zoltan REQUIRED)
#    add_library(Zoltan_Interface INTERFACE)
#    target_link_libraries(Zoltan_Interface INTERFACE zoltan)
#    target_include_directories(Zoltan_Interface INTERFACE ${Zoltan_INCLUDE_DIRS})

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
if (ENABLE_GPU_BUILD)
    set(EXE farrsight_gpu)
    add_executable(${EXE} run_FARRSIGHT.cpp)
    target_link_libraries(${EXE} AMRObj_gpu)
else()
    set(EXE farrsight_cpu)
    add_executable(${EXE} run_FARRSIGHT.cpp)
    target_link_libraries(${EXE} AMRObj_cpu)
endif()

target_link_libraries(${EXE}    PanelObj 
                                icsObj)
target_link_libraries(${EXE} OpenMP::OpenMP_CXX)

install(TARGETS ${EXE} DESTINATION bin)


enable_testing()

add_test(NAME Runs COMMAND ${EXE} 1 0.0 1.0 -2.0 2.0)