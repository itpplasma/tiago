file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "ECE opacity generation failed: ${result}")
endif()
foreach(kernel coordinates coordinates2 o1 o2 o2_blue x1_tenuous x1_flr x1_finite x2 emission)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/ece_opacity_${kernel}.f90"
        "${OUTPUT}/src/generated/ece_opacity_${kernel}.f90" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Regenerated ece_opacity_${kernel}.f90 differs")
    endif()
endforeach()
