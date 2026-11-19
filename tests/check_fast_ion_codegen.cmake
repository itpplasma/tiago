file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "fast-ion kernel generation failed: ${result}")
endif()
foreach(kernel dd_cross_section dd_reactivity cx_hydrogen)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${kernel}.f90"
        "${OUTPUT}/src/generated/${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "regenerated ${kernel} differs from committed source")
    endif()
endforeach()
