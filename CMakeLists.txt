# much of this is owed to Leighton Wilson

cmake_minimum_required(VERSION 3.9...3.14)
project(FARRSIGHT LANGUAGES C CXX)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release" "None")
endif()

set(default_include /usr/local/include)
if(NOT INCLUDE_HEADERS)
    message(STATUS "Setting include directory to '${default_include}' as none was specified.")
    set(INCLUDE_HEADERS ${default_include} CACHE FILEPATH "Choose include directory")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(ENABLE_GPU_BUILD "build OpenACC GPU version"             OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(OpenMP     REQUIRED)
find_package(MPI        REQUIRED)
find_library(BARYTREE_LIB BaryTree_cpu REQUIRED)

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
#include_directories(inc)
#include_directories(inc_eigen/eigen-3.3.7)
#link_directories(src)

# add the executable
add_executable(farrsight_cpu run_FARRSIGHT.cpp)
target_link_libraries     (farrsight_cpu PUBLIC   OpenMP::OpenMP_CXX
                                                MPI::MPI_CXX 
                                                ${BARYTREE_LIB}
                                                FieldObj
                                                PanelObj
                                                icsObj
                                                AMRObj_cpu)
target_include_directories(farrsight_cpu PUBLIC ${INCLUDE_HEADERS})
target_include_directories(farrsight_cpu PUBLIC   inc)

install(TARGETS farrsight_cpu LIBRARY DESTINATION lib
                            ARCHIVE DESTINATION lib)

if (ENABLE_GPU_BUILD)
    find_library(BARYTREE_LIB_GPU BaryTree_gpu REQUIRED)
    target_link_libraries(farrsight_gpu PUBLIC   OpenMP::OpenMP_CXX
                                                MPI::MPI_CXX 
                                                ${BARYTREE_LIB_GPU}
                                                FieldObj_gpu
                                                PanelObj
                                                icsObj
                                                AMRObj_cpu)
    target_include_directories(farrsight_gpu PUBLIC ${INCLUDE_HEADERS})
    target_include_directories(farrsight_gpu PUBLIC   inc)

    install(TARGETS farrsight_gpu LIBRARY DESTINATION lib
                                ARCHIVE DESTINATION lib)
endif()
    
#else()
#    set(EXE farrsight_cpu)
#    add_executable(${EXE} run_FARRSIGHT.cpp)
#    target_link_libraries(${EXE} PUBLIC AMRObj_cpu)
#endif()

#target_link_libraries(${EXE}    PanelObj 
#                                FieldObj
#                                icsObj)
#target_link_libraries(${EXE}    PUBLIC  OpenMP::OpenMP_CXX
#                                        MPI::MPI_CXX
#                                        PanelObj 
#                                        FieldObj
#                                        icsObj)

#install(TARGETS ${EXE} DESTINATION bin)


#enable_testing()

#add_test(NAME Runs COMMAND ${EXE} 1 0.0 1.0 -2.0 2.0)