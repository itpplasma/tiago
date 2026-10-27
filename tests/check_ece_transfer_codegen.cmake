file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "ECE transfer generation failed: ${result}")
endif()
foreach(kernel ece_loss ece_loss_small ece_thermal_cell ece_reflection_loss
        ece_two_stream_wall ece_radiometer_noise)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${kernel}.f90"
        "${OUTPUT}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ${kernel}.f90 differs")
    endif()
endforeach()
