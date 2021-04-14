cmake_minimum_required (VERSION 3.17 FATAL_ERROR)

project (FastVQA CXX)

set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

#
# TO BE SET BY USER
#
if (NOT DEFINED XACC_DIR)
	set(XACC_DIR ~/.xacc)
endif()

message(${XACC_DIR})
	
get_filename_component(DXACC_DIR ${XACC_DIR} ABSOLUTE)
                                              
find_package(XACC REQUIRED)

set (CMAKE_CXX_FLAGS "-g3 -Wall -O2")

add_subdirectory(src)
add_subdirectory(tests)


