cmake_minimum_required(VERSION 3.24)
project(TIAGO LANGUAGES Fortran)

option(TIAGO_ENABLE_NATIVE_OPTIMIZATION "Optimize for the build machine CPU" OFF)
if(TIAGO_ENABLE_NATIVE_OPTIMIZATION)
    if(CMAKE_CROSSCOMPILING)
        message(FATAL_ERROR "Native optimization requires a native build")
    endif()
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
        add_compile_options(-mcpu=native)
    else()
        add_compile_options(-march=native -mtune=native)
    endif()
endif()

include(FetchContent)

set(TIAGO_LIBNEO_GIT "https://github.com/itpplasma/libneo.git" CACHE STRING "libneo upstream to fetch")
set(TIAGO_LIBNEO_TAG "f1bef44e80c4508f04594d4d2ca120eaf17564f8" CACHE STRING "libneo commit/branch/tag used by Tiago")
option(TIAGO_MATCH_LIBNEO_BRANCH "Use the matching libneo branch when available" OFF)
# Use a local effective ref so removing an override restores the pinned default.
set(_tiago_libneo_ref "${TIAGO_LIBNEO_TAG}")
if(DEFINED LIBNEO_BRANCH AND NOT "${LIBNEO_BRANCH}" STREQUAL "")
    set(_tiago_libneo_ref "${LIBNEO_BRANCH}")
elseif(DEFINED ENV{LIBNEO_BRANCH} AND NOT "$ENV{LIBNEO_BRANCH}" STREQUAL "")
    set(_tiago_libneo_ref "$ENV{LIBNEO_BRANCH}")
elseif(TIAGO_MATCH_LIBNEO_BRANCH)
    execute_process(COMMAND git symbolic-ref --quiet --short HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE _tiago_branch OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _tiago_branch_status)
    if(_tiago_branch_status EQUAL 0)
        execute_process(
            COMMAND git ls-remote --heads "${TIAGO_LIBNEO_GIT}" "refs/heads/${_tiago_branch}"
            OUTPUT_VARIABLE _libneo_branch_exists OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _libneo_query_status)
        if(_libneo_query_status EQUAL 0 AND NOT "${_libneo_branch_exists}" STREQUAL "")
            set(_tiago_libneo_ref "${_tiago_branch}")
        endif()
    endif()
endif()
set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo's own tests when vendored" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "Disable libneo OpenMP flags unless explicitly requested" FORCE)
# Tiago uses libneo only to read coil files; its numerics (and its own fortnum) are not needed.
set(LIBNEO_BUILD_NUMERICS OFF CACHE BOOL "Skip libneo numerics components" FORCE)

option(TIAGO_ENABLE_OPENMP "Enable OpenMP parallelism inside Tiago loops" ON)

FetchContent_Declare(libneo
    GIT_REPOSITORY ${TIAGO_LIBNEO_GIT}
    GIT_TAG        ${_tiago_libneo_ref}
    SOURCE_SUBDIR  .)
FetchContent_MakeAvailable(libneo)

add_library(tiago_diagnostics STATIC
    src/diagnostics/diagnostic_types.f90
    src/diagnostics/flux_loops.f90
    src/diagnostics/segmented_rogowski.f90
    src/diagnostics/bprobes.f90)
set_target_properties(tiago_diagnostics PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_diagnostics PUBLIC LIBNEO::field LIBNEO::nctools)

add_library(tiago_solver STATIC
    src/solver/coil_kernels.f90
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

# Lets gfortran vectorise log/sqrt in the Biot-Savart kernels (glibc vector math);
# results are unchanged, only errno is not set.
set_source_files_properties(src/solver/coil_kernels.f90 PROPERTIES COMPILE_OPTIONS -fno-math-errno)

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
    endif()
endif()

add_executable(tiago_cli src/cli/tiago_cli.f90)
target_link_libraries(tiago_cli PRIVATE tiago_diagnostics)

add_library(tiago_build_config STATIC ${TIAGO_GENERATED_DIR}/tiago_config.f90)

add_executable(tiago_vacuum_cli src/cli/tiago_vacuum_cli.f90)
target_link_libraries(tiago_vacuum_cli PRIVATE tiago_solver tiago_build_config)

# Equilibrium reconstruction (tiago_reconstruct): Levenberg-Marquardt in Tiago
# (LAPACK), figures from fortplot, equilibria and their adjoint from VMEC++
# (C++, built below).
option(TIAGO_ENABLE_RECONSTRUCTION
    "Build tiago_reconstruct (fetches fortplot and VMEC++)" OFF)
if(TIAGO_ENABLE_RECONSTRUCTION)
    enable_language(C)
    set(FORTPLOT_BUILD_FORTARRAY_ADAPTER OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(fortplot
        GIT_REPOSITORY https://github.com/lazy-fortran/fortplot.git
        GIT_TAG af5f26e1d35547e9136faec0ba87c1279ae8f04f)
    FetchContent_MakeAvailable(fortplot)
    find_package(LAPACK REQUIRED)

    # VMEC++ (C++) for fixed-boundary equilibria and their exact implicit
    # adjoint. Its exact force Jacobian is generated by Enzyme, a Clang plugin,
    # so the C++ compiler must be Clang with -DTIAGO_ENZYME_PLUGIN pointing at
    # the matching ClangEnzyme-NN.so (see README). VMEC++ installs no headers or
    # CMake package; it is built from a pinned commit, without its Python module
    # (cmake/PatchVmecpp.cmake), and bound to Fortran by src/equilibrium.
    enable_language(CXX)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "tiago_reconstruct builds VMEC++ with Enzyme and needs Clang as "
            "C++ compiler, e.g. -DCMAKE_CXX_COMPILER=clang++-20 "
            "-DTIAGO_ENZYME_PLUGIN=/path/to/ClangEnzyme-20.so, or -DTIAGO_ENABLE_RECONSTRUCTION=OFF")
    endif()
    set(TIAGO_ENZYME_PLUGIN "" CACHE FILEPATH "ClangEnzyme-NN.so matching the Clang C++ compiler")
    if(NOT EXISTS "${TIAGO_ENZYME_PLUGIN}")
        message(FATAL_ERROR "set -DTIAGO_ENZYME_PLUGIN=/path/to/ClangEnzyme-NN.so (see README)")
    endif()
    option(TIAGO_EXPERIMENTAL_ENZYME_TOOLCHAIN
        "Allow an unverified Clang major version for reconstruction adjoints" OFF)
    if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 20 AND
            CMAKE_CXX_COMPILER_VERSION VERSION_LESS 21) AND
       NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 22 AND
            CMAKE_CXX_COMPILER_VERSION VERSION_LESS 23))
        if(NOT TIAGO_EXPERIMENTAL_ENZYME_TOOLCHAIN)
            message(FATAL_ERROR "Reconstruction adjoints are verified with Clang 20.1.8 / "
                "Enzyme v0.0.264 and Clang 22.1.8 / Enzyme v0.0.293. Other Clang majors require "
                "-DTIAGO_EXPERIMENTAL_ENZYME_TOOLCHAIN=ON and independent derivative validation.")
        endif()
        message(WARNING "Experimental reconstruction toolchain: validate the physical adjoints "
            "against re-solved finite differences before using reconstruction results.")
    endif()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -Xclang -load
        -Xclang ${TIAGO_ENZYME_PLUGIN} -x c++ -fsyntax-only /dev/null
        RESULT_VARIABLE _tiago_enzyme_load_result ERROR_VARIABLE _tiago_enzyme_load_error)
    if(NOT _tiago_enzyme_load_result EQUAL 0)
        message(FATAL_ERROR "Enzyme plugin cannot load into the selected Clang compiler: "
            "${_tiago_enzyme_load_error}")
    endif()
    set(TIAGO_VMECPP_HDF5_URL "" CACHE STRING "Mirror of the HDF5 source archive VMEC++ builds")
    set(TIAGO_VMECPP_NETCDF_URL "" CACHE STRING "Mirror of the netCDF-C source archive VMEC++ builds")
    set(VMECPP_ENABLE_ENZYME ON CACHE BOOL "" FORCE)
    set(VMECPP_ENZYME_PLUGIN "${TIAGO_ENZYME_PLUGIN}" CACHE FILEPATH "" FORCE)
    set(VMECPP_HWCAPS_DISPATCH OFF CACHE BOOL "" FORCE)
    # Use LLVM OpenMP for Clang calls and generic/GOMP entry points. Clang has no code
    # generation for libgomp (-fopenmp=libgomp silently drops the pragmas, so
    # VMEC++'s radial threads would not run), while libomp also implements
    # the GOMP ABI of gfortran's OpenMP code. Executables link libomp ahead of
    # libgomp, so its symbols serve both (LLVM also ships libgomp.so.1 as a
    # link to libomp; CMake warns that it hides GCC's libgomp, as intended).
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libomp.so
        OUTPUT_VARIABLE TIAGO_LIBOMP OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT EXISTS "${TIAGO_LIBOMP}")
        message(FATAL_ERROR "libomp for ${CMAKE_CXX_COMPILER} not found (e.g. apt install libomp-20-dev)")
    endif()
    foreach(lang C CXX)
        set(OpenMP_${lang}_FLAGS "-fopenmp=libomp" CACHE STRING "" FORCE)
        set(OpenMP_${lang}_LIB_NAMES "omp" CACHE STRING "" FORCE)
    endforeach()
    set(OpenMP_omp_LIBRARY "${TIAGO_LIBOMP}" CACHE FILEPATH "" FORCE)
    FetchContent_Declare(vmecpp
        GIT_REPOSITORY https://github.com/proximafusion/vmecpp.git
        GIT_TAG 860df56938ea6aaf7a1775c89d51069bab4120b3
        PATCH_COMMAND ${CMAKE_COMMAND} -DSOURCE=<SOURCE_DIR> -DHDF5_URL=${TIAGO_VMECPP_HDF5_URL}
            -DNETCDF_URL=${TIAGO_VMECPP_NETCDF_URL} -P ${CMAKE_SOURCE_DIR}/cmake/PatchVmecpp.cmake
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(vmecpp)
    FetchContent_GetProperties(eigen)
    FetchContent_GetProperties(abseil-cpp)
    FetchContent_GetProperties(abscab-cpp)
    add_library(tiago_vmecpp_adapter STATIC src/equilibrium/vmecpp_adapter.cc)
    target_compile_features(tiago_vmecpp_adapter PRIVATE cxx_std_20)
    target_include_directories(tiago_vmecpp_adapter PRIVATE
        ${vmecpp_SOURCE_DIR}/src/vmecpp/cpp ${eigen_SOURCE_DIR} ${abseil-cpp_SOURCE_DIR}
        ${abscab-cpp_SOURCE_DIR})
    # ABI-relevant: must match VMEC++'s own Eigen alignment (set in its directory scope).
    target_compile_definitions(tiago_vmecpp_adapter PRIVATE
        EIGEN_MAX_ALIGN_BYTES=32 EIGEN_MAX_STATIC_ALIGN_BYTES=32)
    target_link_libraries(tiago_vmecpp_adapter PUBLIC vmecpp_core)
    # VMEC++ reads JSON; classic &INDATA files are converted by its indata2json.
    file(GENERATE OUTPUT ${TIAGO_GENERATED_DIR}/tiago_vmecpp_paths.f90 CONTENT
"module tiago_vmecpp_paths
    implicit none
    character(len=*), parameter :: indata2json_path = &
        \"$<TARGET_FILE:indata2json>\"
end module tiago_vmecpp_paths
")
    set_source_files_properties(${TIAGO_GENERATED_DIR}/tiago_vmecpp_paths.f90 PROPERTIES
        COMPILE_OPTIONS -ffree-line-length-none)
    add_library(tiago_reconstruction STATIC
        ${TIAGO_GENERATED_DIR}/tiago_vmecpp_paths.f90
        src/equilibrium/vmecpp.f90
        src/equilibrium/vmec_edge.f90
        src/equilibrium/equilibrium.f90
        src/reconstruct/reconstruction.f90
        src/reconstruct/reconstruction_plots.f90)
    set_target_properties(tiago_reconstruction PROPERTIES
        Fortran_STANDARD 2018
        Fortran_STANDARD_REQUIRED YES)
    target_link_libraries(tiago_reconstruction PUBLIC ${TIAGO_LIBOMP} tiago_solver
        fortplot::fortplot tiago_vmecpp_adapter ${LAPACK_LIBRARIES})
    # Dependency sorting can place libneo's libgomp before libomp. Force libomp
    # into DT_NEEDED first, so generic omp_* and Clang's __kmpc_* calls use the
    # same runtime. LLVM also supplies the GOMP entry points used by gfortran.
    target_link_options(tiago_reconstruction INTERFACE
        "LINKER:--no-as-needed" "${TIAGO_LIBOMP}" "LINKER:--as-needed")
    set_target_properties(tiago_reconstruction PROPERTIES LINKER_LANGUAGE Fortran)
    add_dependencies(tiago_reconstruction indata2json)
    add_executable(test_vmecpp_equilibrium tests/test_vmecpp_equilibrium.f90)
    target_link_libraries(test_vmecpp_equilibrium PRIVATE tiago_reconstruction)
    add_executable(test_reconstruction_shape_cache tests/test_reconstruction_shape_cache.f90)
    target_link_libraries(test_reconstruction_shape_cache PRIVATE tiago_reconstruction)
    add_executable(tiago_reconstruct src/cli/tiago_reconstruct.f90)
    target_link_libraries(tiago_reconstruct PRIVATE tiago_reconstruction tiago_build_config)
    add_executable(tiago_plot_pulls benchmarks/reconstruction/plot_pulls.f90)
    target_link_libraries(tiago_plot_pulls PRIVATE fortplot::fortplot)
endif()

enable_testing()

add_executable(tiago_coil_loader_tests tests/test_coil_loader.f90)
target_link_libraries(tiago_coil_loader_tests PRIVATE tiago_solver)
add_test(NAME tiago_coil_loader_tests
    COMMAND $<TARGET_FILE:tiago_coil_loader_tests> ${CMAKE_SOURCE_DIR}/tests/data)

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
add_executable(test_plasma_projection tests/test_plasma_projection.f90)
target_link_libraries(test_plasma_projection PRIVATE tiago_solver)
add_test(NAME tiago_plasma_projection COMMAND $<TARGET_FILE:test_plasma_projection>)
add_executable(test_plasma_support tests/test_plasma_support.f90)
target_link_libraries(test_plasma_support PRIVATE tiago_solver)
add_test(NAME tiago_plasma_support COMMAND $<TARGET_FILE:test_plasma_support> ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_plasma_support PROPERTIES FIXTURES_REQUIRED wout_li383)
add_test(NAME tiago_plasma_covariant
    COMMAND $<TARGET_FILE:test_plasma_support> ${TIAGO_WOUT_SIMSOPT_SMALL} covariant)
set_tests_properties(tiago_plasma_covariant PROPERTIES FIXTURES_REQUIRED wout_li383)
add_test(NAME tiago_plasma_conservative
    COMMAND $<TARGET_FILE:test_plasma_support> ${TIAGO_WOUT_SIMSOPT_SMALL} conservative)
set_tests_properties(tiago_plasma_conservative PROPERTIES FIXTURES_REQUIRED wout_li383)



# Vacuum equilibria (no plasma currents): the sheet field must vanish outside.
add_executable(test_plasma_vacuum tests/test_plasma_vacuum.f90)
target_link_libraries(test_plasma_vacuum PRIVATE tiago_solver)
foreach(case IN ITEMS purely_toroidal_field_reference LandremanPaul2021_QA_lowres)
    set(wout ${CMAKE_BINARY_DIR}/tests/cases/wout_${case}.nc)
    add_test(NAME tiago_download_wout_${case}
        COMMAND ${CMAKE_COMMAND}
                -DURL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/wout_${case}.nc
                -DDEST=${wout} -P ${CMAKE_SOURCE_DIR}/cmake/DownloadFile.cmake)
    set_tests_properties(tiago_download_wout_${case} PROPERTIES FIXTURES_SETUP wout_${case} TIMEOUT 120)
    add_test(NAME tiago_plasma_vacuum_${case} COMMAND $<TARGET_FILE:test_plasma_vacuum> ${wout})
    set_tests_properties(tiago_plasma_vacuum_${case} PROPERTIES FIXTURES_REQUIRED wout_${case})
    add_test(NAME tiago_plasma_vacuum_conservative_${case}
        COMMAND $<TARGET_FILE:test_plasma_vacuum> ${wout} conservative)
    set_tests_properties(tiago_plasma_vacuum_conservative_${case} PROPERTIES
        FIXTURES_REQUIRED wout_${case})
endforeach()

# CLI: the README example and argument errors (message, exit status 1)
add_test(NAME tiago_cli_sample
    COMMAND $<TARGET_FILE:tiago_vacuum_cli>
            --coils ${CMAKE_SOURCE_DIR}/tests/data/coil_sample.neo
            --flux ${CMAKE_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --segrog ${CMAKE_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --seg-area 3.40e-4 --samples 8 --output-dir ${CMAKE_BINARY_DIR}/tests/cli_sample)
foreach(bad IN ITEMS
        "no_diagnostics|nothing to evaluate|--coils;x.coils"
        "no_source|no field source|--flux;x.diagno"
        "positional|unknown option|x.coils;x.diagno;y.diagno"
        "bad_samples|--samples expects an integer|--flux;x;--coils;y;--samples;6x"
        "missing_value|--flux requires an argument|--coils;y;--flux")
    string(REPLACE "|" ";" parts "${bad}")
    list(POP_FRONT parts name regex)
    add_test(NAME tiago_cli_rejects_${name} COMMAND $<TARGET_FILE:tiago_vacuum_cli> ${parts})
    set_tests_properties(tiago_cli_rejects_${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "${regex}" WILL_FAIL FALSE)
endforeach()

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

if(TIAGO_ENABLE_RECONSTRUCTION)
    # VMEC++ equilibrium chain (edge field, edge reverse mode, full Jacobian vs
    # re-solved equilibria) and a noise-free reconstruction on LI383 (pinned input).
    set(TIAGO_INPUT_LI383 ${CMAKE_BINARY_DIR}/tests/cases/input.li383_low_res)
    add_test(NAME tiago_download_input_li383
        COMMAND ${CMAKE_COMMAND}
            -DURL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/input.li383_low_res
            -DDEST=${TIAGO_INPUT_LI383} -P ${CMAKE_SOURCE_DIR}/cmake/DownloadFile.cmake)
    set_tests_properties(tiago_download_input_li383 PROPERTIES FIXTURES_SETUP input_li383 TIMEOUT 120)
    add_test(NAME tiago_vmecpp_equilibrium
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${CMAKE_BINARY_DIR}/tests/vmecpp_equilibrium)
    add_test(NAME tiago_vmecpp_parallel
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${CMAKE_BINARY_DIR}/tests/vmecpp_parallel)
    set_tests_properties(tiago_vmecpp_parallel PROPERTIES FIXTURES_REQUIRED input_li383
        TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;TIAGO_VMECPP_VERIFY=1")
    add_test(NAME tiago_reconstruct_exact
        COMMAND ${CMAKE_SOURCE_DIR}/tests/reconstruct_exact.sh $<TARGET_FILE:tiago_reconstruct>
                ${TIAGO_INPUT_LI383} ${CMAKE_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${CMAKE_BINARY_DIR}/tests/reconstruct_exact)
    add_test(NAME tiago_reconstruct_conservative
        COMMAND ${CMAKE_SOURCE_DIR}/tests/reconstruct_exact.sh $<TARGET_FILE:tiago_reconstruct>
                ${TIAGO_INPUT_LI383} ${CMAKE_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${CMAKE_BINARY_DIR}/tests/reconstruct_conservative conservative)
    set_tests_properties(tiago_vmecpp_equilibrium tiago_reconstruct_exact
        tiago_reconstruct_conservative PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600)
    add_test(NAME tiago_reconstruction_shape_cache
        COMMAND $<TARGET_FILE:test_reconstruction_shape_cache> ${TIAGO_INPUT_LI383}
                ${CMAKE_BINARY_DIR}/tests/reconstruction_shape_cache)
    set_tests_properties(tiago_reconstruction_shape_cache PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600)
    add_test(NAME tiago_vmecpp_covariant
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${CMAKE_BINARY_DIR}/tests/vmecpp_covariant 1e-5 covariant)
    set_tests_properties(tiago_vmecpp_covariant PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600 ENVIRONMENT TIAGO_VMECPP_VERIFY=1)
    # check the extracted adjoint blocks against isolated Hessian-vector products
    set_tests_properties(tiago_vmecpp_equilibrium PROPERTIES ENVIRONMENT TIAGO_VMECPP_VERIFY=1)
endif()
