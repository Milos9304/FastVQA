cmake_minimum_required (VERSION 3.15 FATAL_ERROR)

option(COMPILE_TESTS "Compile tests" ON) #ON by default

project (FastVQA
	LANGUAGES C CXX
	DESCRIPTION "High performance variational quantum algorithms simulation framework"
	VERSION 0.1.0
)
set(namespace "fastVQA")

set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

include(GNUInstallDirs)

message(
	STATUS
	"CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
)

set(public_headers
    include/fastVQA.h
    include/vqas/vqe.h
    include/vqas/vqeOptions.h
    include/vqas/vqaOptions.h
    include/vqas/qaoaOptions.h
    include/ansatz.h
    include/cost_function.h
    include/accelerator.h
    include/accelerator_base.h
    include/AqcPqcAccelerator.h
    include/hamiltonian.h
    include/experimentBuffer.h
    include/circuit.h
    include/gate.h
    include/optimizer.h
    include/symbolic_manipulation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/include/QuEST.h
    ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/include/QuEST_precision.h
)
	        
set(quest_src
	quantum_backends/QuEST/QuEST/src/mt19937ar.c
	quantum_backends/QuEST/QuEST/src/QuEST.c
	quantum_backends/QuEST/QuEST/src/QuEST_common.c
	quantum_backends/QuEST/QuEST/src/QuEST_qasm.c
	quantum_backends/QuEST/QuEST/src/QuEST_validation.c
	quantum_backends/QuEST/QuEST/src/CPU/QuEST_cpu.c
	quantum_backends/QuEST/QuEST/src/CPU/QuEST_cpu_local.c
)	        
	                                              
if(DEFINED GMP_H_DIR)
	set (CMAKE_CXX_FLAGS "-I${GMP_H_DIR} -L${GMP_L_DIR} -g3 -Wall -O2")
else()
	set (CMAKE_CXX_FLAGS "-I ~/.local/include -L ~/.local/lib -g3 -Wall -O2")
endif()

add_subdirectory(src/optimizer)

add_library(${PROJECT_NAME} STATIC)
target_sources(${PROJECT_NAME}
    PRIVATE
        src/vqas/vqe.cpp
        src/aqc_pqc/aqc_pqc_accelerator.cpp
        src/gate.cpp
        src/circuit.cpp
        src/ansatz.cpp
        src/accelerator.cpp
        src/hamiltonian.cpp
        src/accelerator_base.cpp
        src/symbolic_manipulation.cpp
        ${quest_src}
        ${optimizer_src}
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/internal
        ${CMAKE_CURRENT_SOURCE_DIR}/include/internal/vqas
        ${CMAKE_CURRENT_SOURCE_DIR}/src/optimizer/nlopt/src/api
        ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/include
        ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/src
        #${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/src/CPU
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")
target_link_libraries(${PROJECT_NAME} PRIVATE nlopt)

install(TARGETS ${PROJECT_NAME} nlopt
   EXPORT "${PROJECT_NAME}Targets"
   PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
   INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# generate and install export file
    install(EXPORT "${PROJECT_NAME}Targets"
        FILE "${PROJECT_NAME}Targets.cmake"
        NAMESPACE ${namespace}::
    DESTINATION cmake
)

include(CMakePackageConfigHelpers)

# generate the version file for the config file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
# create config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION cmake
)
# install config files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION cmake
)
# generate the export targets for the build tree
export(EXPORT "${PROJECT_NAME}Targets"
    FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake"
    NAMESPACE ${namespace}::
)
if(COMPILE_TESTS)
    add_subdirectory(tests)
endif(COMPILE_TESTS)
unset(COMPILE_TESTS CACHE)

