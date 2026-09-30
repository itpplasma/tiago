cmake_minimum_required(VERSION 3.24)
project(TIAGO LANGUAGES Fortran)

# Enable native CPU optimizations
add_compile_options(-march=native -mtune=native)

include(FetchContent)

set(TIAGO_LIBNEO_GIT "git@github.com:itpplasma/libneo.git" CACHE STRING "libneo upstream to fetch")
set(TIAGO_LIBNEO_TAG "f1bef44e80c4508f04594d4d2ca120eaf17564f8" CACHE STRING "libneo commit/branch/tag used by Tiago")
set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo's own tests when vendored" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "Disable libneo OpenMP flags unless explicitly requested" FORCE)

option(TIAGO_ENABLE_OPENMP "Enable OpenMP parallelism inside Tiago loops" ON)

FetchContent_Declare(libneo
    GIT_REPOSITORY ${TIAGO_LIBNEO_GIT}
    GIT_TAG        ${TIAGO_LIBNEO_TAG}
    SOURCE_SUBDIR  .)
FetchContent_MakeAvailable(libneo)

add_library(tiago_core INTERFACE)
target_link_libraries(tiago_core INTERFACE LIBNEO::neo)

add_library(tiago_diagnostics STATIC
    src/diagnostics/diagnostic_types.f90
    src/diagnostics/flux_loops.f90
    src/diagnostics/segmented_rogowski.f90)
set_target_properties(tiago_diagnostics PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_diagnostics PUBLIC tiago_core)

add_library(tiago_solver STATIC
    src/solver/vacuum_forward.f90
    src/solver/coil_loader.f90
    src/solver/plasma_support.f90)
set_target_properties(tiago_solver PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_solver PUBLIC tiago_diagnostics)

# Small VMEC equilibrium fetched on demand by `tiago_vacuum_cli --plasma-sample`
# and used by the plasma regression test.
set(TIAGO_WOUT_SIMSOPT_SMALL_URL
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/wout_li383_low_res_reference.nc")
set(TIAGO_WOUT_SIMSOPT_SMALL ${CMAKE_BINARY_DIR}/tests/cases/wout_li383_low_res_reference.nc)

set(TIAGO_GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${TIAGO_GENERATED_DIR})
configure_file(${CMAKE_SOURCE_DIR}/cmake/tiago_config.f90.in
    ${TIAGO_GENERATED_DIR}/tiago_config.f90 @ONLY)
# Embedded paths/URLs can exceed the 132-column free-form limit.
set_source_files_properties(${TIAGO_GENERATED_DIR}/tiago_config.f90 PROPERTIES
    COMPILE_OPTIONS -ffree-line-length-none)

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

add_executable(tiago_coil_loader_tests tests/test_coil_loader.f90)
target_link_libraries(tiago_coil_loader_tests PRIVATE tiago_solver)
add_test(NAME tiago_coil_loader_tests COMMAND $<TARGET_FILE:tiago_coil_loader_tests>)

add_executable(test_vacuum_semantics tests/test_vacuum_semantics.f90)
target_link_libraries(test_vacuum_semantics PRIVATE tiago_solver)
add_test(NAME tiago_vacuum_semantics
    COMMAND $<TARGET_FILE:test_vacuum_semantics>
            ${CMAKE_SOURCE_DIR}/tests/data/coils_sample.coils)

# Plasma model regression test on a small VMEC equilibrium (downloaded, pinned).
add_test(NAME tiago_download_wout_li383
    COMMAND ${CMAKE_COMMAND} -DURL=${TIAGO_WOUT_SIMSOPT_SMALL_URL} -DDEST=${TIAGO_WOUT_SIMSOPT_SMALL}
            -P ${CMAKE_SOURCE_DIR}/cmake/DownloadFile.cmake)
set_tests_properties(tiago_download_wout_li383 PROPERTIES FIXTURES_SETUP wout_li383 TIMEOUT 120)
add_executable(test_plasma_support tests/test_plasma_support.f90)
target_link_libraries(test_plasma_support PRIVATE tiago_solver)
add_test(NAME tiago_plasma_support COMMAND $<TARGET_FILE:test_plasma_support> ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_plasma_support PROPERTIES FIXTURES_REQUIRED wout_li383)

add_test(NAME tiago_diag_lint_flux
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --kind flux)

add_test(NAME tiago_diag_lint_segrog
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --kind segrog)

add_test(NAME tiago_diag_lint_invalid
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_invalid.diagno
            --kind flux)
set_tests_properties(tiago_diag_lint_invalid PROPERTIES WILL_FAIL TRUE)
