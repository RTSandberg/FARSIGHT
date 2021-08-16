include_directories(${PROJECT_SOURCE_DIR}/inc)


add_library(FieldObj FieldStructure.cpp)
target_link_libraries(FieldObj       PUBLIC  ${BARYTREE_LIB}
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(FieldObj PUBLIC ${INCLUDE_HEADERS})

add_library(icsObj initial_distributions.cpp)
add_library(PanelObj Panel.cpp)

add_library(AMRObj_cpu  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    AMRStructure_field.cpp
                    AMRStructure_step.cpp
                    AMRStructure_interp.cpp)
target_link_libraries(AMRObj_cpu    PUBLIC  PanelObj
                                            icsObj
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
                                                ${INCLUDE_HEADERS})

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

#    set(TRGT AMRObj_gpu)
#    add_library(${TRGT}  AMRStructure.cpp 
#                        AMRStructure_io.cpp
#                        AMRStructure_amr.cpp
#                        AMRStructure_field.cpp
#                        AMRStructure_step.cpp
#                        AMRStructure_interp.cpp)
#    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
#    target_link_libraries(AMRObj_gpu    PUBLIC
#                                        OpenMP::OpenMP_CXX
#                                        MPI::MPI_CXX
#                                        OpenACC_CXX
#                                        $<$<CXX_COMPILER_ID:GNU>:m>)
#    target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
#                                                    ${INCLUDE_HEADERS})
endif()
