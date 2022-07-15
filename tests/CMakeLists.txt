enable_testing()

add_subdirectory(googletest)

include_directories(../include
					../include/vqas
					../include/aqc_pqc
					../quantum_backends/QuEST/QuEST/include/)
                      
add_executable(ansatz_test ansatz_test.cpp)
add_executable(hamiltonian_test hamiltonian_test.cpp)

find_package(FastVQA)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(ansatz_test PRIVATE FastVQA gtest Eigen3::Eigen)
target_link_libraries(hamiltonian_test PRIVATE FastVQA gtest Eigen3::Eigen)

add_test(NAME ansatz_test COMMAND ansatz_test)
add_test(NAME hamiltonian_test COMMAND ansatz_test)
