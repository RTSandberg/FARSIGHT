include_directories(${PROJECT_SOURCE_DIR}/inc)



add_library(SimObj AMRSimulation.cpp
                   AMRSimulation_load.cpp
                   AMRSimulation_io.cpp)
target_link_libraries(SimObj       PUBLIC  AMRObj_cpu
                                            FieldObj
                                            ${BARYTREE_LIB}
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(SimObj PUBLIC ${INCLUDE_HEADERS})