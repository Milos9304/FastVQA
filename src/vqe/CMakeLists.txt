set(LIBRARY_NAME vqe)

file(GLOB SRC vqe.cpp)
set(SRC_ORIG ${SRC})
                            
add_library(${LIBRARY_NAME} SHARED ${SRC})

#add_library(libquest SHARED IMPORTED)
#set_property(TARGET libquest PROPERTY IMPORTED_LOCATION ${DXACC_DIR}/lib/libQuEST.so)

target_link_libraries(${LIBRARY_NAME} PUBLIC io ansatz optimizer)# libquest)#xacc::quantum_gate xacc::quantum_annealing xacc::pauli xacc::pauli xacc::fermion)
