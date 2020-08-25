include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7)
include_directories(/usr/local/include)

set(TRGT AMRObj_cpu)
add_library(PanelObj Panel.cpp)
add_library(${TRGT}  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    AMRStructure_field.cpp
                    AMRStructure_step.cpp
                    AMRStructure_interp.cpp)
add_library(icsObj initial_distributions.cpp)


# target_link_libraries(${TRGT}       PRIVATE BaryTree_lib)

target_link_libraries(${TRGT}       BaryTree_lib
                                    OpenMP::OpenMP_CXX
                                    MPI::MPI_C
                                    $<$<CXX_COMPILER_ID:GNU>:m>)

if(ENABLE_GPU_BUILD)
    set(TRGT AMRObj_gpu)
    add_library(${TRGT}  AMRStructure.cpp 
                        AMRStructure_io.cpp
                        AMRStructure_amr.cpp
                        AMRStructure_field.cpp
                        AMRStructure_step.cpp
                        AMRStructure_interp.cpp)
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_link_libraries(AMRObj_gpu    PUBLIC
                                        OpenMP::OpenMP_CXX
                                        MPI::MPI_C
                                        OpenACC_CXX
                                        $<$<CXX_COMPILER_ID:GNU>:m>)
endif()