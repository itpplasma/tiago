file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Island kernel generation failed: ${result}")
endif()
foreach(kernel fieldline_slope cylindrical_components map_invariants helical_phase)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${OUTPUT}/src/generated/${kernel}.f90"
        "${SOURCE}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ${kernel} differs from committed kernel")
    endif()
endforeach()
