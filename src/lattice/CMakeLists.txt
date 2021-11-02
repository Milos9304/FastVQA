set(LIBRARY_NAME lattice)

file(GLOB SRC *.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})
include_directories(../io)
target_link_libraries(${LIBRARY_NAME} PUBLIC io gmpxx gmp)#xacc::quantum_gate xacc::quantum_annealing xacc::pauli xacc::pauli xacc::fermion)
