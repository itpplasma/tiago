file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "digitiser kernel generation failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SOURCE}/src/generated/digitizer_signal.f90"
    "${OUTPUT}/src/generated/digitizer_signal.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "regenerated digitizer_signal differs from committed source")
endif()
