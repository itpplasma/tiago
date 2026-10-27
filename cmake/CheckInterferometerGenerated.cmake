file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Interferometer generation failed: ${result}")
endif()
foreach(kernel IN ITEMS dispersion_phase dispersion_column dispersion_iq dispersion_atan)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${OUTPUT}/src/generated/${kernel}.f90"
        "${SOURCE}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ${kernel} differs from committed kernel")
    endif()
endforeach()
