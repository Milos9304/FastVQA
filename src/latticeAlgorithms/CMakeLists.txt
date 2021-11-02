set(LIBRARY_NAME latticeAlgorithms)

file(GLOB SRC *.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})