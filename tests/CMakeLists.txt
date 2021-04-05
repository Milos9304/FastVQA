enable_testing()

add_subdirectory(googletest)

add_executable(svpToQubo svpToQubo.cpp ../src/vqaConfig.cpp ../src/symbolic_manipulation.cpp ../src/lattice.cpp)
target_link_libraries(svpToQubo PRIVATE gtest gmpxx gmp)

add_test(NAME svpToQubo COMMAND svpToQubo)