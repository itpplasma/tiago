file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Doppler kernel generation failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SOURCE}/src/generated/doppler_reflectometry.f90"
    "${OUTPUT}/src/generated/doppler_reflectometry.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Generated Doppler kernel differs from the committed source")
endif()
