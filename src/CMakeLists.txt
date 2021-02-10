include_directories(${PROJECT_SOURCE_DIR}/inc)


add_library(FieldObj FieldStructure.cpp)
target_link_libraries(FieldObj       PUBLIC  ${BARYTREE_LIB}
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(FieldObj PUBLIC ${INCLUDE_HEADERS})

if(ENABLE_GPU_BUILD)
    add_library(FieldObj_gpu FieldStructure.cpp)
    target_compile_definitions(FieldObj_gpu PRIVATE OPENACC_ENABLED)
    target_link_libraries(FieldObj_gpu  PUBLIC  ${BARYTREE_LIB_GPU}
                                                OpenMP::OpenMP_CXX
                                                MPI::MPI_CXX
                                                OpenACC_CXX
                                                $<$<CXX_COMPILER_ID:GNU>:m>)
    target_compile_options(FieldObj_gpu PRIVATE -O0)
    target_include_directories(FieldObj_gpu PUBLIC ${INCLUDE_HEADERS})
endif()

add_subdirectory(structure_files)
add_subdirectory(simulation_files)

