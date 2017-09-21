set(ASCGRID AscGrid.cpp)

set(predP_M_FILES predictPoint_mfd.cpp mfdRegionSimilarity.cpp ${ASCGRID})

add_executable(predP_M ${predP_M_FILES})

set(MY_TARGETS predP_M)

foreach (c_target ${MY_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARY})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
endforeach (c_target ${MY_TARGETS})