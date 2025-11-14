cmake_minimum_required(VERSION 3.24)
project(TIAGO LANGUAGES Fortran CXX)

# Enable native CPU optimizations
add_compile_options(-march=native -mtune=native)

include(FetchContent)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(TIAGO_LIBNEO_GIT "git@github.com:itpplasma/libneo.git" CACHE STRING "libneo upstream to fetch")
set(TIAGO_LIBNEO_TAG "main" CACHE STRING "libneo branch/tag used by Tiago")
set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo's own tests when vendored" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "Disable libneo OpenMP flags unless explicitly requested" FORCE)

option(TIAGO_ENABLE_OPENMP "Enable OpenMP parallelism inside Tiago loops" ON)
option(TIAGO_ENABLE_PLASMA "Enable plasma response via virtual casing" ON)

# Find OpenMP early for plasma dependencies
find_package(OpenMP)

FetchContent_Declare(libneo
    GIT_REPOSITORY ${TIAGO_LIBNEO_GIT}
    GIT_TAG        ${TIAGO_LIBNEO_TAG}
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  .)
FetchContent_MakeAvailable(libneo)

if(TIAGO_ENABLE_PLASMA)
    FetchContent_Declare(virtual_casing
        GIT_REPOSITORY "https://github.com/hiddenSymmetries/virtual-casing.git"
        GIT_TAG        "master"
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(virtual_casing)

    # OpenMP is already found above, add to module path for FFTW
    list(APPEND CMAKE_MODULE_PATH ${virtual_casing_SOURCE_DIR}/cmake)

    # Find dependencies for virtual-casing
    find_package(FFTW COMPONENTS DOUBLE_LIB)
    if(NOT FFTW_DOUBLE_LIB_FOUND)
        message(WARNING "FFTW not found, disabling plasma response")
        set(TIAGO_ENABLE_PLASMA OFF)
    endif()
    find_package(BLAS)
    find_package(LAPACK)

    # Create static library from virtual-casing source
    add_library(vc_static STATIC EXCLUDE_FROM_ALL
        ${virtual_casing_SOURCE_DIR}/src/virtual-casing.cpp)
    target_include_directories(vc_static PRIVATE
        ${virtual_casing_SOURCE_DIR}/include
        ${virtual_casing_SOURCE_DIR}/extern/SCTL/include
        ${virtual_casing_SOURCE_DIR}/extern/BIEST/include)
    target_link_libraries(vc_static PRIVATE
        FFTW::Double OpenMP::OpenMP_CXX ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

    set(TIAGO_USE_PLASMA ON)
else()
    set(TIAGO_USE_PLASMA OFF)
endif()

add_library(tiago_core INTERFACE)
target_link_libraries(tiago_core INTERFACE LIBNEO::neo)

add_library(tiago_diagnostics STATIC
    src/diagnostics/diagnostic_types.f90
    src/diagnostics/flux_loops.f90
    src/diagnostics/segmented_rogowski.f90
    src/diagnostics/metadata_registry.f90)
set_target_properties(tiago_diagnostics PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_diagnostics PUBLIC tiago_core)

add_library(tiago_solver STATIC
    src/solver/vacuum_forward.f90
    src/solver/coil_loader.f90)
if(TIAGO_USE_PLASMA)
    target_sources(tiago_solver PRIVATE
        src/solver/plasma_response.f90
        src/solver/plasma_support_enabled.f90)
else()
    target_sources(tiago_solver PRIVATE
        src/solver/plasma_support_stub.f90)
endif()
set_target_properties(tiago_solver PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_solver PUBLIC tiago_diagnostics)
if(TIAGO_USE_PLASMA)
    target_link_libraries(tiago_solver PRIVATE vc_static)
    target_include_directories(tiago_solver PRIVATE
        ${virtual_casing_SOURCE_DIR}/include)
    target_compile_definitions(tiago_solver PUBLIC TIAGO_PLASMA_RESPONSE)
endif()

set(TIAGO_GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${TIAGO_GENERATED_DIR})
configure_file(${CMAKE_SOURCE_DIR}/cmake/tiago_config.f90.in
    ${TIAGO_GENERATED_DIR}/tiago_config.f90 @ONLY)

if(TIAGO_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(tiago_solver PUBLIC OpenMP::OpenMP_Fortran)
        target_compile_definitions(tiago_solver PUBLIC TIAGO_USE_OPENMP)
    endif()
endif()

add_executable(tiago_cli src/cli/tiago_cli.f90)
target_link_libraries(tiago_cli PRIVATE tiago_diagnostics)

add_executable(tiago_vacuum_cli src/cli/tiago_vacuum_cli.f90)
target_sources(tiago_vacuum_cli PRIVATE ${TIAGO_GENERATED_DIR}/tiago_config.f90)
target_link_libraries(tiago_vacuum_cli PRIVATE tiago_solver)

enable_testing()

set(TIAGO_OUTPUT_DIR ${CMAKE_BINARY_DIR}/tests/output)
file(MAKE_DIRECTORY ${TIAGO_OUTPUT_DIR})
set(TIAGO_CASES_DIR ${CMAKE_BINARY_DIR}/tests/cases)
file(MAKE_DIRECTORY ${TIAGO_CASES_DIR})
set(TIAGO_WOUT_LI383_URL "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_li383_low_res_reference.nc")
set(TIAGO_WOUT_LI383 ${TIAGO_CASES_DIR}/wout_li383_low_res_reference.nc)
set(TIAGO_SIMSOPT_DIR ${TIAGO_OUTPUT_DIR}/simsopt_li383)
set(TIAGO_SIMSOPT_GAMMA ${TIAGO_SIMSOPT_DIR}/gamma.csv)
set(TIAGO_SIMSOPT_BTOTAL ${TIAGO_SIMSOPT_DIR}/b_total.csv)
set(TIAGO_SIMSOPT_BEXT ${TIAGO_SIMSOPT_DIR}/b_external.csv)
set(TIAGO_WOUT_SIMSOPT_SMALL_URL
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/tests/test_files/wout_simsopt_small.nc")
set(TIAGO_WOUT_SIMSOPT_SMALL ${CMAKE_SOURCE_DIR}/tests/cases/wout_simsopt_small.nc)
add_executable(tiago_coil_loader_tests tests/test_coil_loader.f90)
target_link_libraries(tiago_coil_loader_tests PRIVATE tiago_solver)
add_test(NAME tiago_coil_loader_tests COMMAND $<TARGET_FILE:tiago_coil_loader_tests>)

if(TIAGO_USE_PLASMA)
    # DISABLED: test_plasma_support_sample requires reference data regeneration
    # add_executable(test_plasma_support_sample tests/test_plasma_support_sample.f90)
    # target_link_libraries(test_plasma_support_sample PRIVATE tiago_solver)
    # add_test(NAME tiago_plasma_support_sample
    #     COMMAND $<TARGET_FILE:test_plasma_support_sample> ${TIAGO_WOUT_LI383})
    # set_tests_properties(tiago_plasma_support_sample PROPERTIES
    #     FIXTURES_REQUIRED wout_li383)
endif()

# DISABLED: test_plasma_response uses point evaluation API not available in virtual-casing
# add_executable(tiago_plasma_response_tests tests/test_plasma_response.f90)
# target_link_libraries(tiago_plasma_response_tests PRIVATE tiago_solver)
# add_test(NAME tiago_plasma_response_tests COMMAND $<TARGET_FILE:tiago_plasma_response_tests>)

add_executable(test_tiago_with_simsopt_grid tests/test_tiago_with_simsopt_grid.f90)
target_link_libraries(test_tiago_with_simsopt_grid PRIVATE tiago_solver)

add_test(NAME tiago_download_wout_li383
    COMMAND ${CMAKE_COMMAND}
            -DURL=${TIAGO_WOUT_LI383_URL}
            -DDEST=${TIAGO_WOUT_LI383}
            -P ${CMAKE_SOURCE_DIR}/cmake/DownloadFile.cmake)
set_tests_properties(tiago_download_wout_li383 PROPERTIES
    FIXTURES_SETUP wout_li383
    TIMEOUT 120)

add_test(NAME tiago_generate_simsopt_li383
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/generate_simsopt_reference.py
            --wout ${TIAGO_WOUT_LI383}
            --output-dir ${TIAGO_SIMSOPT_DIR}
            --src-nphi 16
            --src-ntheta 16)
set_tests_properties(tiago_generate_simsopt_li383 PROPERTIES
    FIXTURES_SETUP simsopt_li383
    FIXTURES_REQUIRED wout_li383)

add_test(NAME tiago_vs_simsopt_reference
    COMMAND $<TARGET_FILE:test_tiago_with_simsopt_grid>
            ${TIAGO_SIMSOPT_GAMMA}
            ${TIAGO_SIMSOPT_BTOTAL}
            ${TIAGO_OUTPUT_DIR}/tiago_b_ext.csv
            ${TIAGO_SIMSOPT_BEXT})
set_tests_properties(tiago_vs_simsopt_reference PROPERTIES
    FIXTURES_REQUIRED simsopt_li383
    FIXTURES_SETUP tiago_bext_fixture)

add_test(NAME tiago_vs_simsopt_plots
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/compare_plasma_vacuum.py
            --tiago-bext ${TIAGO_OUTPUT_DIR}/tiago_b_ext.csv
            --simsopt-btotal ${TIAGO_SIMSOPT_BTOTAL}
            --simsopt-bexternal ${TIAGO_SIMSOPT_BEXT}
            --output-dir ${TIAGO_OUTPUT_DIR}
            --label li383)
set_tests_properties(tiago_vs_simsopt_plots PROPERTIES
    FIXTURES_REQUIRED tiago_bext_fixture)

add_test(NAME tiago_diag_lint_flux
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --kind flux)

add_test(NAME tiago_diag_lint_segrog
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --metadata ${CMAKE_SOURCE_DIR}/docs/diagnostics/registry.json
            --family w7x_segrog_core)

add_test(NAME tiago_diag_lint_invalid
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_invalid.diagno
            --kind flux)
set_tests_properties(tiago_diag_lint_invalid PROPERTIES WILL_FAIL TRUE)


add_test(NAME tiago_vs_xdiagno
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/run_xdiagno.py
            --coil ${CMAKE_SOURCE_DIR}/tests/data/coils_sample.coils
            --flux ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --segrog ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --output ${TIAGO_OUTPUT_DIR}/sample
            --label sample
            --tiago-bin $<TARGET_FILE:tiago_vacuum_cli>
            --seg-area 3.40e-4
            --coil-extcur ${CMAKE_SOURCE_DIR}/tests/data/extcur_sample.input)
set_tests_properties(tiago_vs_xdiagno PROPERTIES SKIP_RETURN_CODE 125)

add_test(NAME tiago_vs_xdiagno_varying_current
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/run_xdiagno.py
            --coil ${CMAKE_SOURCE_DIR}/tests/data/coils_varying.coils
            --flux ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --segrog ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --output ${TIAGO_OUTPUT_DIR}/varying_current
            --label varying_current
            --tiago-bin $<TARGET_FILE:tiago_vacuum_cli>
            --seg-area 3.40e-4
            --coil-extcur ${CMAKE_SOURCE_DIR}/tests/data/extcur_varying.input)
set_tests_properties(tiago_vs_xdiagno_varying_current PROPERTIES SKIP_RETURN_CODE 125)

add_test(NAME tiago_vs_xdiagno_negative_current
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/scripts/run_xdiagno.py
            --coil ${CMAKE_SOURCE_DIR}/tests/data/coils_negative.coils
            --flux ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --segrog ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --output ${TIAGO_OUTPUT_DIR}/negative_current
            --label negative_current
            --tiago-bin $<TARGET_FILE:tiago_vacuum_cli>
            --seg-area 3.40e-4
            --coil-extcur ${CMAKE_SOURCE_DIR}/tests/data/extcur_negative.input)
set_tests_properties(tiago_vs_xdiagno_negative_current PROPERTIES SKIP_RETURN_CODE 125)

foreach(case_name IN ITEMS ncsx_nfp1 ncsx_nfp3)
    if(case_name STREQUAL "ncsx_nfp1")
        set(coil_file "coils.NCSX_nfp1")
        set(coil_url "https://raw.githubusercontent.com/PrincetonUniversity/STELLOPT/develop/BENCHMARKS/FIELDLINES_TEST/coils.NCSX_nfp1")
        set(case_nfp 1)
        set(case_samples 12)
        set(case_tol 5.0e-5)
    else()
        set(coil_file "coils.NCSX")
        set(coil_url "https://raw.githubusercontent.com/PrincetonUniversity/STELLOPT/develop/BENCHMARKS/FIELDLINES_TEST/coils.NCSX")
        set(case_nfp 3)
        set(case_samples 14)
        set(case_tol 5.0e-5)
    endif()
    set(case_rel_tol 1.0e-2)
    add_test(NAME tiago_vs_xdiagno_${case_name}
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/scripts/run_xdiagno.py
                --coil ${CMAKE_SOURCE_DIR}/tests/cases/${case_name}/${coil_file}
                --coil-url ${coil_url}
                --flux ${CMAKE_SOURCE_DIR}/tests/cases/${case_name}/fluxloop.diagno
                --segrog ${CMAKE_SOURCE_DIR}/tests/cases/${case_name}/segrog.diagno
                --output ${TIAGO_OUTPUT_DIR}/${case_name}
                --label ${case_name}
                --tiago-bin $<TARGET_FILE:tiago_vacuum_cli>
                --seg-area 3.40e-4
                --samples ${case_samples}
                --nfp ${case_nfp}
                --tolerance ${case_tol}
                --rel-tolerance ${case_rel_tol}
                --flux-turns ${CMAKE_SOURCE_DIR}/tests/cases/${case_name}/flux_turns.csv
                --segrog-turns ${CMAKE_SOURCE_DIR}/tests/cases/${case_name}/segrog_turns.csv)
    set_tests_properties(tiago_vs_xdiagno_${case_name} PROPERTIES SKIP_RETURN_CODE 125)
endforeach()
