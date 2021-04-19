set(LIBRARY_NAME src)

file(GLOB SRC   lattiQ.cpp 
				vqaConfig.cpp 
				symbolic_manipulation.cpp
				latticeAlgorithms/*.cpp
				utils.cpp
				executionStatistics.cpp
				)
set(SRC_ORIG ${SRC})

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(qaoa)
add_subdirectory(io)
add_subdirectory(lattice)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable(lattiQ ${SRC_ORIG})
#target_link_directories(fastVQA PUBLIC xacc::xacc xacc::quantum_gate)
#target_link_libraries(fastVQA PUBLIC vqe Eigen3::Eigen fplll gmp gmpxx )
target_link_libraries(lattiQ PUBLIC qaoa lattice fplll gmp gmpxx xacc::xacc )#xacc::pauli)