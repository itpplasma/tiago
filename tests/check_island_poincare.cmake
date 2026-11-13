file(MAKE_DIRECTORY "${OUTPUT}")
file(WRITE "${OUTPUT}/coils.wire"
    "periods 1\nbegin filament\nmirror NIL\n0 0 -50 100000\n0 0 50 0 1 wire\nend\n")
file(WRITE "${OUTPUT}/input.nml" "&island
 coil_files = '${OUTPUT}/coils.wire'
 axis_seed = 3.0, 0.0, axis_given = .true.
 o_seed = 4.0, 0.0, map_turns = 1, m = 1, n = 0
 measure_x = .false., measure_width = .false., half_span = 0.1
 poincare_lines = 3, poincare_crossings = 2
 output = '${OUTPUT}/section'
/
")
execute_process(COMMAND "${TIAGO_ISLAND}" "${OUTPUT}/input.nml"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Poincare-only CLI failed: ${stdout} ${stderr}")
endif()
file(READ "${OUTPUT}/section.json" report)
string(JSON ended GET "${report}" end)
string(JSON phase_type TYPE "${report}" axis phase_error)
if(NOT ended OR NOT phase_type STREQUAL "NULL")
    message(FATAL_ERROR "Given-axis JSON must be complete with unknown phase error null")
endif()
file(STRINGS "${OUTPUT}/section-poincare.csv" rows)
list(POP_FRONT rows header)
list(LENGTH rows count)
if(NOT count EQUAL 6)
    message(FATAL_ERROR "Expected six Poincare samples, got ${count}")
endif()
foreach(row IN LISTS rows)
    string(REPLACE "," ";" fields "${row}")
    list(GET fields 0 line)
    list(GET fields 2 r)
    list(GET fields 3 z)
    list(GET fields 4 status)
    if(NOT status EQUAL 0 OR z LESS -1e-12 OR z GREATER 1e-12)
        message(FATAL_ERROR "Purely toroidal wire field must preserve Z: ${row}")
    endif()
    if(line EQUAL 1 AND (r LESS 3.8999999999 OR r GREATER 3.9000000001))
        message(FATAL_ERROR "Wrong first chord point: ${row}")
    elseif(line EQUAL 2 AND (r LESS 3.9999999999 OR r GREATER 4.0000000001))
        message(FATAL_ERROR "Wrong middle chord point: ${row}")
    elseif(line EQUAL 3 AND (r LESS 4.0999999999 OR r GREATER 4.1000000001))
        message(FATAL_ERROR "Wrong last chord point: ${row}")
    endif()
endforeach()

# A toroidal field has neutral periodic points, not an island O point.
file(READ "${OUTPUT}/input.nml" input)
string(REPLACE "measure_width = .false." "measure_width = .true." input "${input}")
file(WRITE "${OUTPUT}/neutral.nml" "${input}")
execute_process(COMMAND "${TIAGO_ISLAND}" "${OUTPUT}/neutral.nml"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(result EQUAL 0)
    message(FATAL_ERROR "Neutral periodic point was incorrectly accepted as an island O point")
endif()
file(READ "${OUTPUT}/section.json" report)
string(JSON status GET "${report}" O status)
if(NOT status EQUAL 4)
    message(FATAL_ERROR "Expected explicit wrong-point-type status, got ${status}")
endif()
