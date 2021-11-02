set(LIBRARY_NAME src)

file(GLOB SRC *.cpp
			  latticeAlgorithms/*.cpp
			  )
			
set(SRC_ORIG ${SRC})

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if (DISTRIBUTED MATCHES 1)
    set(QUACC_SRC_CPU_ARCHITECTURE_DEPENDENT
        ${CMAKE_CURRENT_SOURCE_DIR}/QuaccSharedUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Quacc.cpp
    )
else()
    set(QUACC_SRC_CPU_ARCHITECTURE_DEPENDENT
        ${CMAKE_CURRENT_SOURCE_DIR}/QuaccSharedUtils.cpp
         ${CMAKE_CURRENT_SOURCE_DIR}/Quacc_distributed.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Quacc_distributed_utils.cpp
    )
endif()

set(QUACC_SRC_ARCHITECTURE_DEPENDENT
    ${QUACC_SRC_CPU_ARCHITECTURE_DEPENDENT}
    PARENT_SCOPE
) 

add_subdirectory(QuEST)
add_subdirectory(io)
add_subdirectory(accelerator)
add_subdirectory(optimizer)
add_subdirectory(circuit)
add_subdirectory(ansatz)
add_subdirectory(hamiltonian)
add_subdirectory(qaoa)
add_subdirectory(vqe)
add_subdirectory(lattice)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable(fastVQA ${SRC_ORIG})

target_link_libraries(fastVQA PUBLIC io vqe qaoa hamiltonian accelerator optimizer lattice fplll gmp gmpxx QuEST)
