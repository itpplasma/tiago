# Runs examples/synthetic/li383.nml and checks that every group was written.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(COPY "${EXAMPLE}/li383.nml" "${EXAMPLE}/li383.sensors" "${EXAMPLE}/loops.diagno"
    "${WOUT}" DESTINATION "${WORK}")
foreach(name geometry profiles)
    execute_process(COMMAND "${NCGEN}" -4 -o ${name}.nc "${EXAMPLE}/${name}.cdl"
        WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "ncgen ${name}.cdl failed")
    endif()
endforeach()
execute_process(COMMAND "${SYNTH}" li383.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "tiago_synth failed: ${out}${err}")
endif()
execute_process(COMMAND "${NCGEN}" -version OUTPUT_QUIET ERROR_QUIET)
get_filename_component(bindir "${NCGEN}" DIRECTORY)
find_program(NCDUMP ncdump HINTS "${bindir}")
if(NCDUMP)
    execute_process(COMMAND "${NCDUMP}" -h li383_synthetic.nc WORKING_DIRECTORY "${WORK}"
        OUTPUT_VARIABLE header)
    foreach(group interferometer thomson xics ece mse magnetics)
        if(NOT header MATCHES "group: ${group} ")
            message(FATAL_ERROR "li383_synthetic.nc has no group ${group}")
        endif()
    endforeach()
endif()
