file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Reflectometry kernel generation failed: ${result}")
endif()
foreach(kernel doppler_turning_shift omode_index)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${kernel}.f90"
        "${OUTPUT}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Generated ${kernel} differs from the committed source")
    endif()
endforeach()
