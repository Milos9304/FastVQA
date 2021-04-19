enable_testing()

add_subdirectory(googletest)

include_directories(${DXACC_DIR}/include/xacc/)
include_directories(${DXACC_DIR}/include/quantum/gate)
include_directories(${DXACC_DIR}/include/cppmicroservices4)

add_library(libquest SHARED IMPORTED)
set_property(TARGET libquest PROPERTY IMPORTED_LOCATION ${DXACC_DIR}/lib/libQuEST.so)

add_executable(svpToQubo svpToQubo.cpp 
                         ../src/vqaConfig.cpp 
                         ../src/symbolic_manipulation.cpp
                         )
                         
add_executable(measurementQubitEndianness measurementQubitEndianness.cpp)

add_executable(iterativeLatticeReductionAlgo iterativeLatticeReductionAlgo.cpp 
											 ../src/vqaConfig.cpp 
											 ../src/qaoa/qaoa.cpp
											 ../src/qaoa/qaoaOptions.cpp
											 ../src/executionStatistics.cpp
											 ../src/utils.cpp
											 ../src/latticeAlgorithms/iterativeLatticeReduction.cpp
											 ../src/latticeAlgorithms/latticeAlgorithm.cpp
											 ../src/symbolic_manipulation.cpp)

target_link_libraries(svpToQubo PRIVATE gtest lattice gmpxx gmp fplll)
target_link_libraries(measurementQubitEndianness PRIVATE gtest xacc::xacc libquest)
target_link_libraries(iterativeLatticeReductionAlgo PRIVATE gtest lattice gmpxx gmp xacc::xacc xacc::quantum_gate fplll io libquest)


add_test(NAME svpToQubo COMMAND svpToQubo)
add_test(NAME measurementQubitEndianness COMMAND measurementQubitEndianness)
add_test(NAME iterativeLatticeReductionAlgo COMMAND iterativeLatticeReductionAlgo)
