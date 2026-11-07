file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "CXRS kernel generation failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SOURCE}/src/generated/cxrs_emission.f90"
    "${OUTPUT}/src/generated/cxrs_emission.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "regenerated cxrs_emission differs from committed source")
endif()
