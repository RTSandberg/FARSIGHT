include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7)

add_library(PanelObj Panel.cpp)
add_library(AMRObj AMRStructure.cpp)
add_library(AMRIOObj AMRStructure_io.cpp)
add_library(AMRMeshObj AMRStructure_amr.cpp)
add_library(AMRFieldObj AMRStructure_field.cpp)
add_library(AMRStepObj AMRStructure_step.cpp)
add_library(AMRInterpObj AMRStructure_interp.cpp)



target_link_libraries(AMRFieldObj OpenMP::OpenMP_CXX)
target_link_libraries(AMRInterpObj OpenMP::OpenMP_CXX)