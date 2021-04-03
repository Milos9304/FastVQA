set(LIBRARY_NAME src)

file(GLOB SRC fastVQA.cpp vqaConfig.cpp symbolic_manipulation.cpp lattice.cpp)
set(SRC_ORIG ${SRC})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(qaoa)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable(fastVQA ${SRC_ORIG})
#target_link_directories(fastVQA PUBLIC xacc::xacc xacc::quantum_gate)
target_link_libraries(fastVQA PUBLIC vqe Eigen3::Eigen)
