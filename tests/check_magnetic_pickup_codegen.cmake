file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Magnetic pickup generation failed: ${result}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SOURCE}/src/generated/magnetic_pickup_leaf.f90"
    "${OUTPUT}/src/generated/magnetic_pickup_leaf.f90" RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Regenerated magnetic pickup leaf differs")
endif()
