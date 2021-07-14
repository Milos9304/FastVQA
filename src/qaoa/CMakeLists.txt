set(LIBRARY_NAME qaoa)

file(GLOB SRC qaoa.cpp qaoaOptions.cpp)
set(SRC_ORIG ${SRC})

usfunctiongetresourcesource(TARGET
                            ${LIBRARY_NAME}
                            OUT
                            SRC)
usfunctiongeneratebundleinit(TARGET
                            ${LIBRARY_NAME}
                            OUT
                            SRC)
                            
add_library(${LIBRARY_NAME} SHARED ${SRC})

set(_bundle_name fastVQA)
set_target_properties(${LIBRARY_NAME}
                    PROPERTIES COMPILE_DEFINITIONS
                                US_BUNDLE_NAME=${_bundle_name}
                                US_BUNDLE_NAME
                                ${_bundle_name})
usfunctionembedresources(TARGET
                        ${LIBRARY_NAME}
                        WORKING_DIRECTORY
                        ${CMAKE_CURRENT_SOURCE_DIR}
                        FILES
                        manifest.json)
xacc_configure_plugin_rpath(${LIBRARY_NAME})

#add_library(libquest SHARED IMPORTED)
#set_property(TARGET libquest PROPERTY IMPORTED_LOCATION ${DXACC_DIR}/lib/libQuEST.so)

target_link_libraries(${LIBRARY_NAME} PUBLIC io xacc::xacc xacc::quantum_gate xacc::pauli)# libquest)#xacc::quantum_gate xacc::quantum_annealing xacc::pauli xacc::pauli xacc::fermion)
