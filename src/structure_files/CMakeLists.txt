include_directories(${PROJECT_SOURCE_DIR}/inc)


add_library(icsObj initial_distributions.cpp)
add_library(PanelObj Panel.cpp)

add_library(AMRObj_cpu  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    AMRStructure_interp.cpp)
target_link_libraries(AMRObj_cpu    PUBLIC  PanelObj
                                            icsObj
                                            OpenMP::OpenMP_CXX
                                            MPI::MPI_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
                                                ${INCLUDE_HEADERS})