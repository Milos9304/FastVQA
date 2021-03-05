set(LIBRARY_NAME src)

file(GLOB SRC fastVQA.cpp config_io.cpp)
set(SRC_ORIG ${SRC})

add_subdirectory(qaoa)

#find_package(Boost 1.74 COMPONENTS program_options REQUIRED)
#include_directories( ${Boost_INCLUDE_DIR} )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../../bin")
add_executable(fastVQA ${SRC_ORIG})
#target_link_directories(fastVQA PUBLIC xacc::xacc xacc::quantum_gate)
target_link_libraries(fastVQA PUBLIC vqe)
#target_link_libraries(fastVQA LINK_PUBLIC ${Boost_LIBRARIES} )
