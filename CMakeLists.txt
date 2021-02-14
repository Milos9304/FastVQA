cmake_minimum_required (VERSION 3.17 FATAL_ERROR)

project (FastVQA CXX)


set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

#
# TO BE SET BY USER. LOCAL PATHS NOT SUPPORTED!
#

if (NOT DEFINED ${DXACC_DIR})
	set(DXACC_DIR /home/milos/.xacc)

find_package(XACC REQUIRED)

set (CMAKE_CXX_FLAGS "-g3 -Wall -O2")

add_subdirectory(src)


