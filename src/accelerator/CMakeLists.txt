set(LIBRARY_NAME accelerator)

file(GLOB SRC *.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})

#target_link_libraries(${LIBRARY_NAME} PUBLIC circuit QuEST hamiltonian ansatz)