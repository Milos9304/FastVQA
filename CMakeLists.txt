cmake_minimum_required (VERSION 3.17 FATAL_ERROR)

project (LattiQ CXX)

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

if(DEFINED GMP_H_DIR)
	set (CMAKE_CXX_FLAGS "-I${GMP_H_DIR} -L${GMP_L_DIR} -g3 -Wall -O2")
else()
	set (CMAKE_CXX_FLAGS "-I~/.local/include -g3 -Wall -O2")
endif()

add_subdirectory(src)
add_subdirectory(tests)


