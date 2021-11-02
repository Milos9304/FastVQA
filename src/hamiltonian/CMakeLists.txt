set(LIBRARY_NAME hamiltonian)

file(GLOB SRC hamiltonian.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})

target_link_libraries(${LIBRARY_NAME} PUBLIC circuit)