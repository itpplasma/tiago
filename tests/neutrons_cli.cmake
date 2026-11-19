# Runs tiago_neutrons on a one-cell distribution and a two-row profile and
# checks that the three totals are printed and the beam-target total is positive.
file(MAKE_DIRECTORY "${OUTPUT}")
file(WRITE "${OUTPUT}/dist.txt" "# one cell\n2.0141\n2 1 1\n0.0 1.0\n60000.0 1000.0\n0.5 1.0\n1.0e14 0.0\n")
file(WRITE "${OUTPUT}/profiles.txt" "# s dVds ne Te Ti Zeff\n0.0 13.0 5.0e19 3000.0 3000.0 1.5\n1.0 13.0 2.0e19 300.0 300.0 1.5\n")
execute_process(COMMAND "${CLI}" "${OUTPUT}/dist.txt" "${OUTPUT}/profiles.txt" 7
    RESULT_VARIABLE result OUTPUT_VARIABLE out)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "tiago_neutrons failed: ${result}")
endif()
foreach(key thermonuclear beam_target beam_beam)
    if(NOT out MATCHES "total_${key}_per_s")
        message(FATAL_ERROR "missing total_${key}_per_s in\n${out}")
    endif()
endforeach()
if(NOT out MATCHES "total_beam_target_per_s +[1-9]")
    message(FATAL_ERROR "beam-target total not positive in\n${out}")
endif()
