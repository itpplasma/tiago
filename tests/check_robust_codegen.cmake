file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Robust signal generation failed: ${result}")
endif()
foreach(kernel student_residual student_residual_small)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${kernel}.f90"
        "${OUTPUT}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ${kernel}.f90 differs")
    endif()
endforeach()
