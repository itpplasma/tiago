execute_process(COMMAND "${CMAKE_COMMAND}"
    -S "${SOURCE}/tests/subproject" -B "${OUTPUT}" -G "${GENERATOR}"
    "-DTIAGO_SOURCE=${SOURCE}"
    "-DCMAKE_C_COMPILER=${C_COMPILER}" "-DCMAKE_Fortran_COMPILER=${FORTRAN_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
    "-DFETCHCONTENT_SOURCE_DIR_FORTNUM=${FORTNUM_SOURCE}"
    "-DFETCHCONTENT_SOURCE_DIR_LIBNEO=${LIBNEO_SOURCE}"
    "-DFETCHCONTENT_SOURCE_DIR_FORTIO=${FORTIO_SOURCE}"
    "-DFETCHCONTENT_SOURCE_DIR_FORTGEN_RUNTIME=${FORTGEN_SOURCE}"
    "-DFETCHCONTENT_SOURCE_DIR_JSONFORTRAN=${JSONFORTRAN_SOURCE}"
    RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "TIAGO parent configure failed: ${status}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build "${OUTPUT}"
    --config "${CONFIG}" --target tiago_consumer --parallel 4 RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "TIAGO parent build failed: ${status}")
endif()
execute_process(COMMAND "${CMAKE_CTEST_COMMAND}" --test-dir "${OUTPUT}"
    -C "${CONFIG}" -R "^consumer$" --output-on-failure RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "TIAGO parent consumer failed: ${status}")
endif()
