enable_testing()

add_subdirectory(googletest)

include_directories(${DXACC_DIR}/include/xacc/)

add_library(libquest SHARED IMPORTED)
set_property(TARGET libquest PROPERTY IMPORTED_LOCATION ${DXACC_DIR}/lib/libQuEST.so)

add_executable(svpToQubo svpToQubo.cpp ../src/vqaConfig.cpp ../src/symbolic_manipulation.cpp ../src/lattice.cpp)
add_executable(measurementQubitEndianness measurementQubitEndianness.cpp)

target_link_libraries(svpToQubo PRIVATE gtest gmpxx gmp fplll)
target_link_libraries(measurementQubitEndianness PRIVATE gtest xacc::xacc libquest)

add_test(NAME svpToQubo COMMAND svpToQubo)
add_test(NAME measurementQubitEndianness COMMAND measurementQubitEndianness)

