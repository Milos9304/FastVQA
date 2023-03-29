cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

project (QuEST
	LANGUAGES C #CXX
	DESCRIPTION "QuEST"
	VERSION 0.1.0
)
set(namespace "quest")

#set(CMAKE_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

set(CMAKE_BUILD_TYPE Release)


#set(quest_headers
#    ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/include/QuEST.h
#    ${CMAKE_CURRENT_SOURCE_DIR}/quantum_backends/QuEST/QuEST/include/QuEST_precision.h
#)
	        
set(quest_src
	QuEST/QuEST/src/mt19937ar.c
	QuEST/QuEST/src/QuEST.c
	QuEST/QuEST/src/QuEST_common.c
	QuEST/QuEST/src/QuEST_qasm.c
	QuEST/QuEST/src/QuEST_validation.c
	QuEST/QuEST/src/CPU/QuEST_cpu.c
	QuEST/QuEST/src/CPU/QuEST_cpu_local.c
)

add_library(${PROJECT_NAME} STATIC)
target_sources(${PROJECT_NAME}
    PRIVATE
        ${quest_src}
)
target_include_directories(${PROJECT_NAME}
    PRIVATE
        QuEST/QuEST/include
        QuEST/QuEST/src
)
target_compile_definitions(${PROJECT_NAME}
    PRIVATE
    QuEST_PREC=${QUEST_PRECISION}
    CMAKE_BUILD
)
