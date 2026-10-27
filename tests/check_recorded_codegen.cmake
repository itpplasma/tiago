file(MAKE_DIRECTORY "${WORK}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "FortSym generator failed: ${status}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${EXPECTED}" "${WORK}/src/generated/recorded_affine.f90"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "recorded_affine differs from the committed FortSym output")
endif()
