file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "ECE recorded generation failed: ${result}")
endif()
foreach(kernel ece_launch_scale ece_field_norm ece_segment_length ece_parallel_fraction
        ece_doppler_x1 ece_doppler_x2)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${kernel}.f90"
        "${OUTPUT}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ${kernel}.f90 differs")
    endif()
endforeach()
