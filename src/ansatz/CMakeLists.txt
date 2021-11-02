set(LIBRARY_NAME ansatz)

file(GLOB SRC ansatz.cpp)
set(SRC_ORIG ${SRC})

add_library(${LIBRARY_NAME} SHARED ${SRC})

target_link_libraries(${LIBRARY_NAME} PUBLIC circuit)