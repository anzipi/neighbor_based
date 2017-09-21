set(ASCGRID AscGrid.cpp)
set(UTILS common_func.cpp UtilsClass/utils.cpp UtilsClass/ModelException.cpp)
INCLUDE_DIRECTORIES(UtilsClass)

set(predP_S_FILES predictPoint_sfd.cpp sfdRegionSimilarity.cpp ${ASCGRID} ${UTILS})
set(predP_M_FILES predictPoint_mfd.cpp mfdRegionSimilarity.cpp ${ASCGRID} ${UTILS})

add_executable(predP_S ${predP_S_FILES})
add_executable(predP_M ${predP_M_FILES})

set(MY_TARGETS predP_S
               predP_M)

foreach (c_target ${MY_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARY})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
endforeach (c_target ${MY_TARGETS})
