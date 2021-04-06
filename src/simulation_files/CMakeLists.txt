include_directories(${PROJECT_SOURCE_DIR}/inc)



add_library(SimObj AMRSimulation.cpp
                   AMRSimulation_load.cpp
                   AMRSimulation_gather_scatter.cpp
                   AMRSimulation_field.cpp
                   AMRSimulation_step.cpp
                   AMRSimulation_remesh.cpp
                   AMRSimulation_run.cpp
                   AMRSimulation_io.cpp)
target_link_libraries(SimObj       PUBLIC  AMRObj_cpu
                                            FieldObj
                                            ${BARYTREE_LIB}
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(SimObj PUBLIC ${INCLUDE_HEADERS})

if (ENABLE_GPU_BUILD)

    add_library(SimObj_gpu AMRSimulation.cpp
                    AMRSimulation_load.cpp
                    AMRSimulation_gather_scatter.cpp
                    AMRSimulation_field.cpp
                    AMRSimulation_step.cpp
                    AMRSimulation_remesh.cpp
                    AMRSimulation_run.cpp
                    AMRSimulation_io.cpp)
    target_compile_definitions(SimObj_gpu PRIVATE OPENACC_ENABLED)
    target_link_libraries(SimObj_gpu       PUBLIC  AMRObj_cpu
                                                FieldObj_gpu
                                                ${BARYTREE_LIB}
                                                OpenMP::OpenMP_CXX
                                                MPI::MPI_CXX
                                                OpenACC_CXX
                                                $<$<CXX_COMPILER_ID:GNU>:m>)
    target_include_directories(SimObj_gpu PUBLIC ${INCLUDE_HEADERS})

endif()