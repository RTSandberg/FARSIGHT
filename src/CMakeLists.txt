include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7)

add_library(PanelObj Panel.cpp)
add_library(AMRObj  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    AMRStructure_field.cpp
                    AMRStructure_step.cpp
                    AMRStructure_interp.cpp)
add_library(icsObj initial_distributions.cpp)


target_link_libraries(AMRObj OpenMP::OpenMP_CXX)