# In the future, ASCGRID will be totally replaced by COMMON_TAUDEM.
# The TauDEM source code is copy from https://github.com/lreis2415/TauDEM
set(ASCGRID AscGrid.cpp)
set(COMMON_TAUDEM TauDEM/commonLib.cpp TauDEM/tiffIO.cpp)
set(UTILS common_func.cpp UtilsClass/utils.cpp UtilsClass/ModelException.cpp)
INCLUDE_DIRECTORIES(UtilsClass TauDEM)

set(predP_S_FILES predictPoint_sfd.cpp sfdRegionSimilarity.cpp ${ASCGRID} ${UTILS} ${COMMON_TAUDEM})
set(predP_M_FILES predictPoint_mfd.cpp mfdRegionSimilarity.cpp ${ASCGRID} ${UTILS} ${COMMON_TAUDEM})

add_executable(predP_S ${predP_S_FILES})
add_executable(predP_M ${predP_M_FILES})

set(MY_TARGETS predP_S
               predP_M)

# TauDEM programs
set(D8 TauDEM/D8FlowDirmn.cpp TauDEM/d8.cpp TauDEM/Node.cpp ${COMMON_TAUDEM})
set(DINF TauDEM/DinfFlowDirmn.cpp TauDEM/dinf.cpp TauDEM/Node.cpp ${COMMON_TAUDEM})
set(PITREMOVE TauDEM/flood.cpp TauDEM/PitRemovemn.cpp ${COMMON_TAUDEM})

add_executable(d8flowdir ${D8})
add_executable(dinfflowdir ${DINF})
add_executable(pitremove ${PITREMOVE})

set(TAUDEM_TARGETS d8flowdir
        dinfflowdir
        pitremove
        )
# TauDEM programs
foreach (c_target ${TAUDEM_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARY})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "TauDEM")
    ENDIF ()
endforeach (c_target ${TAUDEM_TARGETS})
# third-party programs based on TauDEM framework
foreach (c_target ${MY_TARGETS})
    target_link_libraries(${c_target} ${MPI_LIBRARIES} ${GDAL_LIBRARY})
    install(TARGETS ${c_target} DESTINATION ${INSTALL_DIR})
    IF (MSVC)
        SET_PROPERTY(TARGET ${c_target} PROPERTY FOLDER "Neighbor_based")
    ENDIF ()
endforeach (c_target ${MY_TARGETS})
