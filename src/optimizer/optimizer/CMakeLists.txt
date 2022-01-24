set(LIBRARY_NAME optimizer)

file(GLOB SRC *.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})

SET(NLOPT_PYTHON OFF CACHE BOOL "")
SET(NLOPT_OCTAVE OFF CACHE BOOL "")
SET(NLOPT_MATLAB OFF CACHE BOOL "")
SET(NLOPT_GUILE OFF CACHE BOOL "")
SET(NLOPT_SWIG OFF CACHE BOOL "")

add_subdirectory(nlopt)
target_include_directories(
  ${LIBRARY_NAME}
  PUBLIC . .. ${CMAKE_BUILD_DIR}/src/optimizer/nlopt/src/api)

target_link_libraries(${LIBRARY_NAME} PUBLIC nlopt)