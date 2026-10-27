file(MAKE_DIRECTORY "${OUTPUT}/src/generated")
execute_process(COMMAND "${GENERATOR}" WORKING_DIRECTORY "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Thomson kernel regeneration failed: ${result}")
endif()
foreach(name ts_depolarization ts_spectral ts_bessel_chain ts_photoelectrons
             ts_recorded ts_rayleigh ts_noise)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${SOURCE}/src/generated/${name}.f90" "${OUTPUT}/src/generated/${name}.f90"
        RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Generated Thomson kernel differs: ${name}")
    endif()
endforeach()
