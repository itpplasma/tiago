file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "wall generator failed: ${status}")
endif()
foreach(name wall_legacy_energy wall_energy wall_flux wall_thermal wall_radiance)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${name}.f90"
        "${OUTPUT}/src/generated/${name}.f90" RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "wall kernel differs: ${name}")
    endif()
endforeach()
