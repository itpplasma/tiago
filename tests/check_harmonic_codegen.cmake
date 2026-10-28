file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Harmonic mixer generation failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${OUTPUT}/src/generated/harmonic_mixer.f90"
    "${SOURCE}/src/generated/harmonic_mixer.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Regenerated harmonic mixer differs from committed kernel")
endif()
