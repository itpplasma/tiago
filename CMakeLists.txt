cmake_minimum_required(VERSION 3.24)
project(TIAGO LANGUAGES C Fortran)

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
set(TIAGO_LIBNEO_TAG "42b04688a2fa2805bb84273c68dc35e78b50f72e" CACHE STRING "libneo commit/branch/tag used by Tiago")
option(TIAGO_ENABLE_SPECTRE "Build the libneo SPECTRE plasma source (HDF5, fortnum)" OFF)
option(TIAGO_MATCH_LIBNEO_BRANCH "Use the matching libneo branch when available" OFF)
# Use a local effective ref so removing an override restores the pinned default.
set(_tiago_libneo_ref "${TIAGO_LIBNEO_TAG}")
if(DEFINED LIBNEO_BRANCH AND NOT "${LIBNEO_BRANCH}" STREQUAL "")
    set(_tiago_libneo_ref "${LIBNEO_BRANCH}")
elseif(DEFINED ENV{LIBNEO_BRANCH} AND NOT "$ENV{LIBNEO_BRANCH}" STREQUAL "")
    set(_tiago_libneo_ref "$ENV{LIBNEO_BRANCH}")
elseif(TIAGO_MATCH_LIBNEO_BRANCH)
    execute_process(COMMAND git symbolic-ref --quiet --short HEAD
        WORKING_DIRECTORY "${TIAGO_SOURCE_DIR}"
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
# The SPECTRE field and coordinates are part of libneo's numerics targets.
set(LIBNEO_BUILD_NUMERICS ${TIAGO_ENABLE_SPECTRE} CACHE BOOL "Build libneo field evaluators" FORCE)

option(TIAGO_ENABLE_OPENMP "Enable OpenMP parallelism inside Tiago loops" ON)

# fortnum: B-spline profile bases, NUTS and PSIS for the inference core. Declared
# before libneo, whose numerics reuse an existing fortnum target.
set(TIAGO_FORTNUM_TAG "ef26a189d39078200f981cbec62bed2515df493a" CACHE STRING "fortnum commit used by Tiago")
if(NOT TARGET fortnum)
    FetchContent_Declare(fortnum
        GIT_REPOSITORY https://github.com/lazy-fortran/fortnum.git
        GIT_TAG        ${TIAGO_FORTNUM_TAG})
    FetchContent_MakeAvailable(fortnum)
endif()

FetchContent_Declare(libneo
    GIT_REPOSITORY ${TIAGO_LIBNEO_GIT}
    GIT_TAG        ${_tiago_libneo_ref}
    SOURCE_SUBDIR  .)
FetchContent_MakeAvailable(libneo)

add_library(tiago_diagnostics STATIC
    src/diagnostics/diagnostic_types.f90
    src/diagnostics/flux_loops.f90
    src/diagnostics/segmented_rogowski.f90
    src/diagnostics/bprobes.f90
    src/generated/mse_rows.f90
    src/generated/mse_local.f90
    src/diagnostics/mse.f90
    src/diagnostics/sensor_set.f90
    src/generated/hermite_basis.f90
    src/generated/circular_state.f90
    src/generated/cold_resonance.f90
    src/generated/hermite_higher.f90
    src/diagnostics/radial_interp.f90
    src/diagnostics/real_space_state.f90
    src/diagnostics/analytic_state.f90)
set_target_properties(tiago_diagnostics PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_diagnostics PUBLIC LIBNEO::field LIBNEO::nctools)

add_library(tiago_solver STATIC
    src/solver/coil_kernels.f90
    src/solver/vacuum_forward.f90
    src/solver/coil_loader.f90
    src/solver/plasma_field_source.f90
    src/solver/plasma_support.f90
    src/solver/sensor_signals.f90
    src/generated/vmec_r_x.f90
    src/generated/vmec_r_y.f90
    src/generated/vmec_z.f90
    src/generated/vmec_cos.f90
    src/generated/vmec_vector_component.f90
    src/solver/vmec_state.f90)
if(TIAGO_ENABLE_SPECTRE)
    if(NOT TARGET LIBNEO::spectre_field)
        message(FATAL_ERROR "SPECTRE requires libneo's LIBNEO::spectre_field target")
    endif()
    target_sources(tiago_solver PRIVATE src/solver/spectre_source.f90)
    target_link_libraries(tiago_solver PUBLIC LIBNEO::spectre_field)
endif()
set_target_properties(tiago_solver PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_solver PUBLIC tiago_diagnostics)

# Small VMEC equilibrium fetched on demand by `tiago_vacuum_cli --plasma-sample`
# and used by the plasma regression test.
set(TIAGO_WOUT_SIMSOPT_SMALL_URL
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/wout_li383_low_res_reference.nc")
set(TIAGO_WOUT_SIMSOPT_SMALL ${TIAGO_BINARY_DIR}/tests/cases/wout_li383_low_res_reference.nc)

# Lets gfortran vectorise log/sqrt in the Biot-Savart kernels (glibc vector math);
# results are unchanged, only errno is not set.
set_source_files_properties(src/solver/coil_kernels.f90 PROPERTIES COMPILE_OPTIONS -fno-math-errno)

set(TIAGO_GENERATED_DIR ${TIAGO_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${TIAGO_GENERATED_DIR})
configure_file(${TIAGO_SOURCE_DIR}/cmake/tiago_config.f90.in
    ${TIAGO_GENERATED_DIR}/tiago_config.f90 @ONLY)
# Embedded paths/URLs can exceed the 132-column free-form limit.
set_source_files_properties(${TIAGO_GENERATED_DIR}/tiago_config.f90 PROPERTIES
    COMPILE_OPTIONS -ffree-line-length-none)

if(TIAGO_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(tiago_solver PUBLIC OpenMP::OpenMP_Fortran)
        target_link_libraries(tiago_diagnostics PUBLIC OpenMP::OpenMP_Fortran)
    endif()
endif()

# Recorded windows use the standard netCDF C hyperslab API through ISO C.
find_library(TIAGO_NETCDF_C_LIBRARY NAMES netcdf REQUIRED)
find_package(CURL 7.66 REQUIRED)
find_package(Threads REQUIRED)
set(JSONFORTRAN_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(JSONFORTRAN_STATIC_LIBRARY_ONLY ON CACHE BOOL "" FORCE)
set(SKIP_DOC_GEN ON CACHE BOOL "" FORCE)
FetchContent_Declare(jsonfortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG b7dc6b3227e16f69ec1061c7aa9eae4f1d2de943)
FetchContent_MakeAvailable(jsonfortran)
add_library(tiago_signals STATIC
    src/generated/recorded_affine.f90
    src/generated/harmonic_mixer.f90
    src/signals/harmonic_lockin.f90
    src/signals/http_curl.c
    src/signals/http_client.f90
    src/signals/rest_adapter.f90
    src/signals/online_source.f90
    src/signals/signal_source.f90
    src/signals/netcdf_source.f90
    src/generated/integrator_signal.f90
    src/generated/digitizer_signal.f90
    src/signals/adc.f90
    src/signals/integrated_drift.f90
    src/signals/signal_reduction.f90
    src/signals/calibration.f90
    src/signals/sha256.f90
    src/signals/netcdf_c.f90)
set_target_properties(tiago_signals PROPERTIES
    Fortran_STANDARD 2018 Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_signals PUBLIC LIBNEO::nctools fortnum
    PRIVATE ${TIAGO_NETCDF_C_LIBRARY} CURL::libcurl Threads::Threads jsonfortran-static)
if(TIAGO_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(tiago_signals PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Inference core without a provider: kinetic profiles and their prior, the
# profile-diagnostic contract and the reference point probe.
add_library(tiago_inference STATIC
    src/generated/kinetic_pressure.f90
    src/generated/student_residual.f90
    src/generated/student_residual_small.f90
    src/inference/robust_signals.f90
    src/inference/calibration_campaign.f90
    src/inference/kinetic_profiles.f90
    src/inference/profile_diagnostic.f90
    src/inference/point_probe.f90
    src/inference/recorded_diagnostic.f90
    src/inference/integrator_diagnostic.f90
    src/inference/digitizer_diagnostic.f90
    src/inference/magnetic_pickup.f90
    src/generated/magnetic_pickup_leaf.f90
    src/inference/magnetic_response.f90
    src/generated/diamagnetic_energy.f90
    src/generated/thermal_energy_density.f90
    src/inference/energy_consistency.f90
    src/generated/xics_emissivity.f90
    src/inference/sight_geometry.f90
    src/generated/medium_pullback.f90
    src/generated/density_label_jet.f90
    src/inference/state_sampler.f90
    src/generated/vmec_product_pair.f90
    src/inference/vmecpp_state.f90
    src/inference/thomson_scattering.f90
    src/inference/profile_geometry_io.f90
    src/inference/mse_state.f90
    src/inference/point_table.f90
    src/generated/ts_depolarization.f90
    src/generated/ts_spectral.f90
    src/generated/ts_bessel_chain.f90
    src/generated/ts_photoelectrons.f90
    src/generated/ts_recorded.f90
    src/generated/ts_rayleigh.f90
    src/generated/ts_noise.f90
    src/inference/thomson_spectrum.f90
    src/inference/thomson_polychromator.f90
    src/signals/thomson_reduction.f90
    src/generated/ece_hamiltonian_o.f90
    src/generated/ece_hamiltonian_x.f90
    src/generated/ece_hamiltonian_unmagnetized.f90
    src/inference/ece_ray.f90
    src/inference/ece_state_medium.f90
    src/generated/ece_loss.f90
    src/generated/ece_loss_small.f90
    src/generated/ece_thermal_cell.f90
    src/generated/ece_reflection_loss.f90
    src/generated/ece_two_stream_wall.f90
    src/generated/ece_radiometer_noise.f90
    src/inference/ece_transfer.f90
    src/generated/ece_opacity_coordinates.f90
    src/generated/ece_opacity_coordinates2.f90
    src/generated/ece_opacity_emission.f90
    src/generated/ece_opacity_o1.f90
    src/generated/ece_opacity_o2.f90
    src/generated/ece_opacity_o2_blue.f90
    src/generated/ece_opacity_x1_finite.f90
    src/generated/ece_opacity_x1_flr.f90
    src/generated/ece_opacity_x1_tenuous.f90
    src/generated/ece_opacity_x2.f90
    src/inference/ece_opacity.f90
    src/inference/ece_opacity_transfer.f90
    src/inference/ece_receiver.f90
    src/generated/ece_launch_scale.f90
    src/generated/ece_field_norm.f90
    src/generated/ece_segment_length.f90
    src/generated/ece_parallel_fraction.f90
    src/generated/ece_doppler_x1.f90
    src/generated/ece_doppler_x2.f90
    src/inference/ece_radiometer.f90
    src/inference/ece_recorded.f90
    src/inference/line_integrals.f90
    src/generated/dispersion_phase.f90
    src/generated/dispersion_column.f90
    src/generated/dispersion_iq.f90
    src/generated/dispersion_atan.f90
    src/inference/dispersion_interferometer.f90
    src/signals/iq_phase_reduction.f90)
set_target_properties(tiago_inference PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_inference PUBLIC fortnum tiago_diagnostics tiago_signals tiago_solver)
if(TIAGO_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(tiago_inference PUBLIC OpenMP::OpenMP_Fortran)
endif()

add_executable(tiago_cli src/cli/tiago_cli.f90)
target_link_libraries(tiago_cli PRIVATE tiago_diagnostics)

add_library(tiago_build_config STATIC ${TIAGO_GENERATED_DIR}/tiago_config.f90)

# Source revision for output provenance, regenerated on every build.
add_custom_target(tiago_revision_file
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_GENERATED_DIR}/tiago_revision.f90
        -P ${TIAGO_SOURCE_DIR}/cmake/TiagoRevision.cmake
    BYPRODUCTS ${TIAGO_GENERATED_DIR}/tiago_revision.f90)
add_library(tiago_synthetic STATIC ${TIAGO_GENERATED_DIR}/tiago_revision.f90
    src/inference/synthetic_forward.f90)
add_dependencies(tiago_synthetic tiago_revision_file)
set_target_properties(tiago_synthetic PROPERTIES
    Fortran_STANDARD 2018 Fortran_STANDARD_REQUIRED YES)
target_link_libraries(tiago_synthetic PUBLIC tiago_inference tiago_solver
    ${TIAGO_NETCDF_C_LIBRARY})
add_executable(tiago_synth src/cli/tiago_synth.f90)
target_link_libraries(tiago_synth PRIVATE tiago_synthetic)

add_executable(tiago_vacuum_cli src/cli/tiago_vacuum_cli.f90)
set_source_files_properties(src/cli/tiago_vacuum_cli.f90 PROPERTIES Fortran_PREPROCESS ON)
if(TIAGO_ENABLE_SPECTRE)
    target_compile_definitions(tiago_vacuum_cli PRIVATE TIAGO_ENABLE_SPECTRE)
endif()
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
    if(APPLE)
        set(_tiago_omp_name libomp.dylib)
    else()
        set(_tiago_omp_name libomp.so)
    endif()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=${_tiago_omp_name}
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
            -DNETCDF_URL=${TIAGO_VMECPP_NETCDF_URL} -P ${TIAGO_SOURCE_DIR}/cmake/PatchVmecpp.cmake
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(vmecpp)
    # VMEC++ builds Release with ThinLTO. Enzyme >= 0.0.303 with Clang 23 then
    # defers its pass to link time and leaves __enzyme_* unresolved unless lld
    # loads LLDEnzyme; compiling the two Enzyme translation units without LTO
    # runs the pass at compile time with every toolchain.
    set(_tiago_enzyme_arch "")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        set(_tiago_enzyme_arch "-march=x86-64")
    endif()
    set_source_files_properties(
        ${vmecpp_SOURCE_DIR}/src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/exact_force_jvp.cc
        ${vmecpp_SOURCE_DIR}/src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/exact_force_vjp.cc
        DIRECTORY ${vmecpp_SOURCE_DIR}
        PROPERTIES COMPILE_OPTIONS "${_tiago_enzyme_arch};-O2;-fno-lto;-fplugin=${TIAGO_ENZYME_PLUGIN}")
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
        src/reconstruct/reconstruction_plots.f90
        src/inference/joint_model.f90
        src/inference/inference_engines.f90
        src/inference/synthetic_case.f90)
    set_target_properties(tiago_reconstruction PROPERTIES
        Fortran_STANDARD 2018
        Fortran_STANDARD_REQUIRED YES)
    target_link_libraries(tiago_reconstruction PUBLIC ${TIAGO_LIBOMP} tiago_solver tiago_inference
        fortplot::fortplot tiago_vmecpp_adapter ${LAPACK_LIBRARIES})
    # Dependency sorting can place libneo's libgomp before libomp. Force libomp
    # into DT_NEEDED first, so generic omp_* and Clang's __kmpc_* calls use the
    # same runtime. LLVM also supplies the GOMP entry points used by gfortran.
    if(APPLE)
        target_link_options(tiago_reconstruction INTERFACE "${TIAGO_LIBOMP}")
    else()
        target_link_options(tiago_reconstruction INTERFACE
            "LINKER:--no-as-needed" "${TIAGO_LIBOMP}" "LINKER:--as-needed")
    endif()
    set_target_properties(tiago_reconstruction PROPERTIES LINKER_LANGUAGE Fortran)
    add_dependencies(tiago_reconstruction indata2json)
    add_executable(test_vmecpp_equilibrium tests/test_vmecpp_equilibrium.f90)
    target_link_libraries(test_vmecpp_equilibrium PRIVATE tiago_reconstruction)
    add_executable(test_joint_model tests/test_joint_model.f90)
    target_link_libraries(test_joint_model PRIVATE tiago_reconstruction)
    add_executable(test_reconstruction_shape_cache tests/test_reconstruction_shape_cache.f90)
    target_link_libraries(test_reconstruction_shape_cache PRIVATE tiago_reconstruction)
    add_executable(tiago_infer src/cli/tiago_infer.f90)
    target_link_libraries(tiago_infer PRIVATE tiago_reconstruction)
    add_executable(tiago_reconstruct src/cli/tiago_reconstruct.f90)
    target_link_libraries(tiago_reconstruct PRIVATE tiago_reconstruction tiago_build_config)
    add_executable(tiago_plot_pulls benchmarks/reconstruction/plot_pulls.f90)
    target_link_libraries(tiago_plot_pulls PRIVATE fortplot::fortplot)
endif()

# Export module search paths to parent projects linking the core libraries.
set(TIAGO_MODULE_DIR "${TIAGO_BINARY_DIR}/modules")
file(MAKE_DIRECTORY "${TIAGO_MODULE_DIR}")
foreach(target IN ITEMS tiago_diagnostics tiago_solver tiago_signals tiago_inference)
    set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY "${TIAGO_MODULE_DIR}")
    target_include_directories(${target} PUBLIC "$<BUILD_INTERFACE:${TIAGO_MODULE_DIR}>")
endforeach()

enable_testing()

if(TIAGO_ENABLE_SPECTRE)
    add_executable(test_spectre_source tests/test_spectre_source.f90)
    target_link_libraries(test_spectre_source PRIVATE tiago_solver)
    set(spectre_fixture ${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5)
    set(spectre_reference ${libneo_SOURCE_DIR}/test/spectre/data/reference_eval.dat)
    add_test(NAME tiago_spectre_source COMMAND $<TARGET_FILE:test_spectre_source>
        ${spectre_fixture} ${spectre_reference})
    foreach(mode IN ITEMS reject_exterior reject_interior reject_missing_flux)
        add_test(NAME tiago_spectre_${mode} COMMAND $<TARGET_FILE:test_spectre_source>
            ${spectre_fixture} ${spectre_reference} ${mode})
    endforeach()
    set_tests_properties(tiago_spectre_reject_exterior PROPERTIES
        PASS_REGULAR_EXPRESSION "virtual casing requires exterior diagnostic points")
    set_tests_properties(tiago_spectre_reject_interior PROPERTIES
        PASS_REGULAR_EXPRESSION "cannot locate interior SPECTRE point")
    set_tests_properties(tiago_spectre_reject_missing_flux PROPERTIES
        PASS_REGULAR_EXPRESSION "idia=1 requires signed linked plasma flux metadata")
    add_test(NAME tiago_cli_spectre COMMAND $<TARGET_FILE:tiago_vacuum_cli>
        --plasma-spectre ${spectre_fixture}
        --bprobes ${TIAGO_SOURCE_DIR}/tests/data/spectre/probes.diagno
        --plasma-nphi 32 --plasma-ntheta 32
        --output-dir ${TIAGO_BINARY_DIR}/tests/cli_spectre)
    set(nested_wout ${TIAGO_BINARY_DIR}/tests/cases/wout_qa_nested.nc)
    add_test(NAME tiago_download_wout_qa_nested COMMAND ${CMAKE_COMMAND}
        -DURL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc
        -DDEST=${nested_wout} -P ${TIAGO_SOURCE_DIR}/cmake/DownloadFile.cmake)
    set_tests_properties(tiago_download_wout_qa_nested PROPERTIES
        FIXTURES_SETUP qa_nested TIMEOUT 120)
    add_executable(test_spectre_vmec tests/test_spectre_vmec.f90)
    target_link_libraries(test_spectre_vmec PRIVATE tiago_solver)
    add_test(NAME tiago_spectre_vmec COMMAND $<TARGET_FILE:test_spectre_vmec>
        ${TIAGO_SOURCE_DIR}/tests/data/spectre/qa_nested.h5 ${nested_wout}
        ${TIAGO_SOURCE_DIR}/tests/data/spectre/qa_vmec_cartesian.dat)
    add_test(NAME tiago_spectre_vmec_refined COMMAND $<TARGET_FILE:test_spectre_vmec>
        ${TIAGO_SOURCE_DIR}/tests/data/spectre/qa_nested.h5 ${nested_wout}
        ${TIAGO_SOURCE_DIR}/tests/data/spectre/qa_vmec_cartesian.dat 128)
    set_tests_properties(tiago_spectre_vmec tiago_spectre_vmec_refined PROPERTIES
        FIXTURES_REQUIRED qa_nested TIMEOUT 180)
else()
    add_test(NAME tiago_cli_spectre_disabled COMMAND $<TARGET_FILE:tiago_vacuum_cli>
        --plasma-spectre absent.h5 --bprobes absent.diagno)
    set_tests_properties(tiago_cli_spectre_disabled PROPERTIES
        PASS_REGULAR_EXPRESSION "requires a build with TIAGO_ENABLE_SPECTRE=ON")
endif()

add_executable(test_mse tests/test_mse.f90)
target_link_libraries(test_mse PRIVATE tiago_solver)
add_test(NAME tiago_mse COMMAND $<TARGET_FILE:test_mse>
    ${TIAGO_SOURCE_DIR}/tests/data ${TIAGO_SOURCE_DIR}/examples/mse_channels.txt)

add_executable(tiago_coil_loader_tests tests/test_coil_loader.f90)
target_link_libraries(tiago_coil_loader_tests PRIVATE tiago_solver)
add_test(NAME tiago_coil_loader_tests
    COMMAND $<TARGET_FILE:tiago_coil_loader_tests> ${TIAGO_SOURCE_DIR}/tests/data)

add_executable(test_inactive_coil_segments tests/test_inactive_coil_segments.f90)
target_link_libraries(test_inactive_coil_segments PRIVATE tiago_solver)
add_test(NAME tiago_inactive_coil_segments COMMAND $<TARGET_FILE:test_inactive_coil_segments>)

add_executable(test_coil_potential tests/test_coil_potential.f90)
target_link_libraries(test_coil_potential PRIVATE tiago_solver)
add_test(NAME tiago_coil_potential COMMAND $<TARGET_FILE:test_coil_potential>)

add_executable(test_vacuum_semantics tests/test_vacuum_semantics.f90)
target_link_libraries(test_vacuum_semantics PRIVATE tiago_solver)
add_test(NAME tiago_vacuum_semantics
    COMMAND $<TARGET_FILE:test_vacuum_semantics>
            ${TIAGO_SOURCE_DIR}/tests/data/coils_sample.coils)

# Plasma model regression test on a small VMEC equilibrium (downloaded, pinned).
add_test(NAME tiago_download_wout_li383
    COMMAND ${CMAKE_COMMAND} -DURL=${TIAGO_WOUT_SIMSOPT_SMALL_URL} -DDEST=${TIAGO_WOUT_SIMSOPT_SMALL}
            -P ${TIAGO_SOURCE_DIR}/cmake/DownloadFile.cmake)
set_tests_properties(tiago_download_wout_li383 PROPERTIES FIXTURES_SETUP wout_li383 TIMEOUT 120)
add_executable(test_mse_plasma tests/test_mse_plasma.f90)
target_link_libraries(test_mse_plasma PRIVATE tiago_solver)
foreach(mode IN ITEMS contravariant covariant conservative)
    add_test(NAME tiago_mse_plasma_${mode} COMMAND $<TARGET_FILE:test_mse_plasma>
        ${TIAGO_WOUT_SIMSOPT_SMALL} ${TIAGO_SOURCE_DIR}/tests/data/mse_coils.coils ${mode})
    set_tests_properties(tiago_mse_plasma_${mode} PROPERTIES FIXTURES_REQUIRED wout_li383)
endforeach()

add_executable(test_mse_state tests/test_mse_state.f90)
target_link_libraries(test_mse_state PRIVATE tiago_inference)
add_test(NAME tiago_mse_state COMMAND $<TARGET_FILE:test_mse_state> ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_mse_state PROPERTIES FIXTURES_REQUIRED wout_li383)

add_executable(test_sensor_set tests/test_sensor_set.f90)
target_link_libraries(test_sensor_set PRIVATE tiago_solver)
set(TIAGO_TEST_SENSOR_GEOMETRY ${TIAGO_BINARY_DIR}/tests/sensors)
file(MAKE_DIRECTORY ${TIAGO_TEST_SENSOR_GEOMETRY})
add_test(NAME tiago_sensor_set COMMAND $<TARGET_FILE:test_sensor_set>
    ${TIAGO_SOURCE_DIR}/tests/data/sensors/synthetic.sensors)
add_test(NAME tiago_sensor_set_plasma COMMAND $<TARGET_FILE:test_sensor_set>
    ${TIAGO_SOURCE_DIR}/tests/data/sensors/synthetic.sensors ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_sensor_set tiago_sensor_set_plasma PROPERTIES
    ENVIRONMENT TIAGO_TEST_SENSOR_GEOMETRY=${TIAGO_TEST_SENSOR_GEOMETRY} RESOURCE_LOCK sensor_geometry)
set_tests_properties(tiago_sensor_set_plasma PROPERTIES FIXTURES_REQUIRED wout_li383)

add_executable(test_profile_instruments tests/test_profile_instruments.f90)
target_link_libraries(test_profile_instruments PRIVATE tiago_inference tiago_solver)
add_test(NAME tiago_profile_instruments COMMAND $<TARGET_FILE:test_profile_instruments>)
add_test(NAME tiago_profile_instruments_vmec
    COMMAND $<TARGET_FILE:test_profile_instruments> ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_profile_instruments_vmec PROPERTIES FIXTURES_REQUIRED wout_li383)

add_executable(test_synthetic_forward tests/test_synthetic_forward.f90)
target_link_libraries(test_synthetic_forward PRIVATE tiago_synthetic)
set(TIAGO_SYNTH_TEST_DIR ${TIAGO_BINARY_DIR}/tests/synthetic)
file(MAKE_DIRECTORY ${TIAGO_SYNTH_TEST_DIR})
add_test(NAME tiago_synthetic_forward
    COMMAND $<TARGET_FILE:test_synthetic_forward> ${TIAGO_SYNTH_TEST_DIR})
add_test(NAME tiago_synthetic_forward_vmec
    COMMAND $<TARGET_FILE:test_synthetic_forward> ${TIAGO_SYNTH_TEST_DIR} ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_synthetic_forward_vmec PROPERTIES FIXTURES_REQUIRED wout_li383)
add_executable(test_point_table tests/test_point_table.f90)
target_link_libraries(test_point_table PRIVATE tiago_synthetic)
add_test(NAME tiago_point_table
    COMMAND $<TARGET_FILE:test_point_table> ${TIAGO_SYNTH_TEST_DIR} ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_point_table PROPERTIES FIXTURES_REQUIRED wout_li383)
# The public example: CDL inputs -> netCDF-4 (ncgen) -> tiago_synth.
find_program(TIAGO_NCGEN ncgen)
if(TIAGO_NCGEN)
    add_test(NAME tiago_synth_example
        COMMAND ${CMAKE_COMMAND} -DNCGEN=${TIAGO_NCGEN} -DSYNTH=$<TARGET_FILE:tiago_synth>
            -DEXAMPLE=${TIAGO_SOURCE_DIR}/examples/synthetic -DWOUT=${TIAGO_WOUT_SIMSOPT_SMALL}
            -DWORK=${TIAGO_BINARY_DIR}/tests/synthetic_example
            -P ${TIAGO_SOURCE_DIR}/tests/synth_example.cmake)
    set_tests_properties(tiago_synth_example PROPERTIES FIXTURES_REQUIRED wout_li383)
endif()

add_executable(test_plasma_projection tests/test_plasma_projection.f90)
target_link_libraries(test_plasma_projection PRIVATE tiago_solver)
add_test(NAME tiago_plasma_projection COMMAND $<TARGET_FILE:test_plasma_projection>)
add_executable(test_plasma_support tests/test_plasma_support.f90)
target_link_libraries(test_plasma_support PRIVATE tiago_solver)
add_test(NAME tiago_plasma_support COMMAND $<TARGET_FILE:test_plasma_support> ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_plasma_support PROPERTIES FIXTURES_REQUIRED wout_li383)
add_executable(test_field_source tests/test_field_source.f90)
target_link_libraries(test_field_source PRIVATE tiago_solver)
add_test(NAME tiago_field_source COMMAND $<TARGET_FILE:test_field_source>
    ${TIAGO_SOURCE_DIR}/tests/data/coils_sample.coils ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_field_source PROPERTIES FIXTURES_REQUIRED wout_li383)
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
    set(wout ${TIAGO_BINARY_DIR}/tests/cases/wout_${case}.nc)
    add_test(NAME tiago_download_wout_${case}
        COMMAND ${CMAKE_COMMAND}
                -DURL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/wout_${case}.nc
                -DDEST=${wout} -P ${TIAGO_SOURCE_DIR}/cmake/DownloadFile.cmake)
    set_tests_properties(tiago_download_wout_${case} PROPERTIES FIXTURES_SETUP wout_${case} TIMEOUT 120)
    add_test(NAME tiago_plasma_vacuum_${case} COMMAND $<TARGET_FILE:test_plasma_vacuum> ${wout})
    set_tests_properties(tiago_plasma_vacuum_${case} PROPERTIES FIXTURES_REQUIRED wout_${case})
    add_test(NAME tiago_plasma_vacuum_conservative_${case}
        COMMAND $<TARGET_FILE:test_plasma_vacuum> ${wout} conservative)
    set_tests_properties(tiago_plasma_vacuum_conservative_${case} PROPERTIES
        FIXTURES_REQUIRED wout_${case})
endforeach()

add_test(NAME tiago_mse_cli
    COMMAND sh ${TIAGO_SOURCE_DIR}/tests/mse_cli.sh $<TARGET_FILE:tiago_vacuum_cli>
        ${TIAGO_SOURCE_DIR}/tests/data/mse_coils.coils
        ${TIAGO_SOURCE_DIR}/examples/mse_channels.txt ${TIAGO_BINARY_DIR}/tests/mse_cli)

# CLI: the README example and argument errors (message, exit status 1)
add_test(NAME tiago_cli_sample
    COMMAND $<TARGET_FILE:tiago_vacuum_cli>
            --coils ${TIAGO_SOURCE_DIR}/tests/data/coil_sample.neo
            --flux ${TIAGO_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --segrog ${TIAGO_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --seg-area 3.40e-4 --samples 8 --output-dir ${TIAGO_BINARY_DIR}/tests/cli_sample)
foreach(bad IN ITEMS
        "no_diagnostics|nothing to evaluate|--coils;x.coils"
        "no_source|no field source|--flux;x.diagno"
        "positional|unknown option|x.coils;x.diagno;y.diagno"
        "bad_samples|--samples expects an integer|--flux;x;--coils;y;--samples;6x"
        "missing_value|--flux requires an argument|--coils;y;--flux"
        "conflicting_plasma|mutually exclusive|--flux;x;--plasma-wout;y;--plasma-spectre;z"
        "linked_flux_without_source|--plasma-linked-flux requires|--flux;x;--coils;y;--plasma-linked-flux;1")
    string(REPLACE "|" ";" parts "${bad}")
    list(POP_FRONT parts name regex)
    add_test(NAME tiago_cli_rejects_${name} COMMAND $<TARGET_FILE:tiago_vacuum_cli> ${parts})
    set_tests_properties(tiago_cli_rejects_${name} PROPERTIES
        PASS_REGULAR_EXPRESSION "${regex}" WILL_FAIL FALSE)
endforeach()

add_executable(benchmark_recorded EXCLUDE_FROM_ALL
    benchmarks/signals/benchmark_recorded.f90)
target_link_libraries(benchmark_recorded PRIVATE tiago_inference)

add_executable(test_signal_sources tests/test_signal_sources.f90
    tests/grouped_signal_fixture.c)
target_link_libraries(test_signal_sources PRIVATE tiago_signals)
add_test(NAME tiago_signal_sources COMMAND $<TARGET_FILE:test_signal_sources>)
add_executable(test_recorded_diagnostic tests/test_recorded_diagnostic.f90)
target_link_libraries(test_recorded_diagnostic PRIVATE tiago_inference)
add_test(NAME tiago_recorded_diagnostic COMMAND $<TARGET_FILE:test_recorded_diagnostic>)

add_executable(test_profile_diagnostic tests/test_profile_diagnostic.f90)
target_link_libraries(test_profile_diagnostic PRIVATE tiago_inference)
add_test(NAME tiago_profile_diagnostic COMMAND $<TARGET_FILE:test_profile_diagnostic>)
add_executable(test_kinetic_profiles tests/test_kinetic_profiles.f90)
target_link_libraries(test_kinetic_profiles PRIVATE tiago_inference)
add_test(NAME tiago_kinetic_profiles COMMAND $<TARGET_FILE:test_kinetic_profiles>)

add_test(NAME tiago_diag_lint_flux
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${TIAGO_SOURCE_DIR}/tests/data/fluxloop_sample.diagno
            --kind flux)

add_test(NAME tiago_diag_lint_segrog
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${TIAGO_SOURCE_DIR}/tests/data/segrog_sample.diagno
            --kind segrog)

add_test(NAME tiago_diag_lint_invalid
    COMMAND $<TARGET_FILE:tiago_cli> diag lint
            ${TIAGO_SOURCE_DIR}/tests/data/fluxloop_invalid.diagno
            --kind flux)
set_tests_properties(tiago_diag_lint_invalid PROPERTIES WILL_FAIL TRUE)

# Private machine geometry (e.g. W7-X in plasma/data) is never in this
# repository; set TIAGO_W7X_SENSORS to its sensor-set manifest to lint it.
add_test(NAME tiago_w7x_sensor_set
    COMMAND sh -c "test -n \"$TIAGO_W7X_SENSORS\" || exit 77; exec \"$0\" diag lint \"$TIAGO_W7X_SENSORS\" --kind sensors"
        $<TARGET_FILE:tiago_cli>)
set_tests_properties(tiago_w7x_sensor_set PROPERTIES SKIP_RETURN_CODE 77)

if(TIAGO_ENABLE_RECONSTRUCTION)
    # VMEC++ equilibrium chain (edge field, edge reverse mode, full Jacobian vs
    # re-solved equilibria) and a noise-free reconstruction on LI383 (pinned input).
    set(TIAGO_INPUT_LI383 ${TIAGO_BINARY_DIR}/tests/cases/input.li383_low_res)
    add_test(NAME tiago_download_input_li383
        COMMAND ${CMAKE_COMMAND}
            -DURL=https://raw.githubusercontent.com/hiddenSymmetries/simsopt/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/input.li383_low_res
            -DDEST=${TIAGO_INPUT_LI383} -P ${TIAGO_SOURCE_DIR}/cmake/DownloadFile.cmake)
    set_tests_properties(tiago_download_input_li383 PROPERTIES FIXTURES_SETUP input_li383 TIMEOUT 120)
    add_test(NAME tiago_vmecpp_equilibrium
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${TIAGO_BINARY_DIR}/tests/vmecpp_equilibrium)
    add_test(NAME tiago_vmecpp_parallel
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${TIAGO_BINARY_DIR}/tests/vmecpp_parallel)
    set_tests_properties(tiago_vmecpp_parallel PROPERTIES FIXTURES_REQUIRED input_li383
        TIMEOUT 600 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;TIAGO_VMECPP_VERIFY=1")
    add_test(NAME tiago_reconstruct_exact
        COMMAND ${TIAGO_SOURCE_DIR}/tests/reconstruct_exact.sh $<TARGET_FILE:tiago_reconstruct>
                ${TIAGO_INPUT_LI383} ${TIAGO_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${TIAGO_BINARY_DIR}/tests/reconstruct_exact)
    add_test(NAME tiago_reconstruct_conservative
        COMMAND ${TIAGO_SOURCE_DIR}/tests/reconstruct_exact.sh $<TARGET_FILE:tiago_reconstruct>
                ${TIAGO_INPUT_LI383} ${TIAGO_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${TIAGO_BINARY_DIR}/tests/reconstruct_conservative conservative)
    set_tests_properties(tiago_vmecpp_equilibrium tiago_reconstruct_exact
        tiago_reconstruct_conservative PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600)
    add_test(NAME tiago_joint_model
        COMMAND $<TARGET_FILE:test_joint_model> ${TIAGO_INPUT_LI383}
                ${TIAGO_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${TIAGO_BINARY_DIR}/tests/joint_model)
    set_tests_properties(tiago_joint_model PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    add_test(NAME tiago_infer_smoke
        COMMAND ${TIAGO_SOURCE_DIR}/tests/infer_smoke.sh $<TARGET_FILE:tiago_infer> ${TIAGO_INPUT_LI383}
                ${TIAGO_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
                ${TIAGO_BINARY_DIR}/tests/infer_smoke)
    set_tests_properties(tiago_infer_smoke PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 900 ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1")
    add_test(NAME tiago_reconstruction_shape_cache
        COMMAND $<TARGET_FILE:test_reconstruction_shape_cache> ${TIAGO_INPUT_LI383}
                ${TIAGO_BINARY_DIR}/tests/reconstruction_shape_cache)
    set_tests_properties(tiago_reconstruction_shape_cache PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600)
    add_test(NAME tiago_vmecpp_covariant
        COMMAND $<TARGET_FILE:test_vmecpp_equilibrium> ${TIAGO_INPUT_LI383}
                ${TIAGO_BINARY_DIR}/tests/vmecpp_covariant 1e-5 covariant)
    set_tests_properties(tiago_vmecpp_covariant PROPERTIES
        FIXTURES_REQUIRED input_li383 TIMEOUT 600 ENVIRONMENT TIAGO_VMECPP_VERIFY=1)
    # check the extracted adjoint blocks against isolated Hessian-vector products
    set_tests_properties(tiago_vmecpp_equilibrium PROPERTIES ENVIRONMENT TIAGO_VMECPP_VERIFY=1)
endif()

add_executable(test_thomson tests/test_thomson.f90)
target_link_libraries(test_thomson PRIVATE tiago_inference)
add_test(NAME tiago_thomson COMMAND $<TARGET_FILE:test_thomson>)
add_executable(test_recorded_magnetics tests/test_recorded_magnetics.f90)
target_link_libraries(test_recorded_magnetics PRIVATE tiago_inference)
add_test(NAME tiago_recorded_magnetics COMMAND $<TARGET_FILE:test_recorded_magnetics>
    ${TIAGO_SOURCE_DIR}/tests/data/mse_coils.coils ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_recorded_magnetics PROPERTIES FIXTURES_REQUIRED wout_li383)
add_executable(test_integrated_drift tests/test_integrated_drift.f90)
target_link_libraries(test_integrated_drift PRIVATE tiago_signals)
add_test(NAME tiago_integrated_drift COMMAND $<TARGET_FILE:test_integrated_drift>)
add_executable(test_digitizer tests/test_digitizer.f90)
target_link_libraries(test_digitizer PRIVATE tiago_inference)
add_test(NAME tiago_digitizer COMMAND $<TARGET_FILE:test_digitizer>)

set(TIAGO_FORTSYM_BUILD "" CACHE PATH "FortSym build for kernel regeneration checks")
set(FORTSYM_BUILD "" CACHE PATH "FortSym build for optional kernel regeneration tests")
if(TIAGO_FORTSYM_BUILD)
    set(FORTSYM_BUILD "${TIAGO_FORTSYM_BUILD}")
endif()
if(FORTSYM_BUILD)
    add_subdirectory(gen codegen)
    add_test(NAME tiago_ece_transfer_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ece_transfer> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/ece-transfer-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_ece_transfer_codegen.cmake)
    add_test(NAME tiago_thomson_regeneration COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_thomson>
        -DSOURCE=${TIAGO_SOURCE_DIR} -DOUTPUT=${TIAGO_BINARY_DIR}/thomson-regeneration
        -P ${TIAGO_SOURCE_DIR}/tests/thomson_regeneration.cmake)
    add_test(NAME tiago_magnetic_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_magnetics> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/magnetic-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_magnetic_codegen.cmake)
    add_test(NAME tiago_digitizer_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_digitizer> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/digitizer-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_digitizer_codegen.cmake)
    add_test(NAME tiago_interferometer_generated COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_interferometer>
        -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/interferometer_codegen
        -P ${TIAGO_SOURCE_DIR}/cmake/CheckInterferometerGenerated.cmake)
    add_test(NAME tiago_recorded_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_recorded>
        -DEXPECTED=${TIAGO_SOURCE_DIR}/src/generated/recorded_affine.f90
        -DWORK=${TIAGO_BINARY_DIR}/tests/codegen -P
        ${TIAGO_SOURCE_DIR}/tests/check_recorded_codegen.cmake)
endif()
if(TIAGO_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    add_executable(benchmark_magnetics EXCLUDE_FROM_ALL
        benchmarks/signals/benchmark_magnetics.f90)
    target_link_libraries(benchmark_magnetics PRIVATE tiago_inference)
endif()

add_executable(test_ece_ray tests/test_ece_ray.f90 tests/ece_ray_fixture.f90
    src/generated/ece_test_medium.f90)
target_link_libraries(test_ece_ray PRIVATE tiago_inference)
add_test(NAME tiago_ece_ray COMMAND $<TARGET_FILE:test_ece_ray>)

add_executable(benchmark_ece_ray EXCLUDE_FROM_ALL benchmarks/ece/benchmark_ray.f90
    tests/ece_ray_fixture.f90 src/generated/ece_test_medium.f90)
target_link_libraries(benchmark_ece_ray PRIVATE tiago_inference)

add_test(NAME tiago_ece_codegen COMMAND sh ${TIAGO_SOURCE_DIR}/tests/ece_codegen.sh
    ${TIAGO_SOURCE_DIR} "${TIAGO_FORTSYM_BUILD}" ${CMAKE_Fortran_COMPILER})
set_tests_properties(tiago_ece_codegen PROPERTIES SKIP_RETURN_CODE 77 LABELS codegen)

add_executable(test_dispersion_interferometer tests/test_dispersion_interferometer.f90)
target_link_libraries(test_dispersion_interferometer PRIVATE tiago_inference)
add_test(NAME tiago_dispersion_interferometer COMMAND test_dispersion_interferometer)
add_executable(test_interferometer_coverage tests/test_interferometer_coverage.f90)
target_link_libraries(test_interferometer_coverage PRIVATE tiago_inference)
add_test(NAME tiago_interferometer_coverage COMMAND test_interferometer_coverage)

if(TIAGO_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    add_executable(benchmark_interferometer EXCLUDE_FROM_ALL
        benchmarks/interferometer/benchmark_interferometer.f90)
    target_link_libraries(benchmark_interferometer PRIVATE tiago_inference)
endif()
add_executable(test_interferometer_profiles tests/test_interferometer_profiles.f90)
target_link_libraries(test_interferometer_profiles PRIVATE tiago_inference)
add_test(NAME tiago_interferometer_profiles COMMAND test_interferometer_profiles)
# Synthetic loopback server; production code has no Python dependency.
find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(TIAGO_MOCK_DIR ${TIAGO_BINARY_DIR}/tests/mock_archive)
add_test(NAME tiago_mock_start COMMAND ${Python3_EXECUTABLE}
    ${TIAGO_SOURCE_DIR}/tests/mock_archive.py start ${TIAGO_MOCK_DIR})
set_tests_properties(tiago_mock_start PROPERTIES FIXTURES_SETUP archive_http)
add_test(NAME tiago_mock_stop COMMAND ${Python3_EXECUTABLE}
    ${TIAGO_SOURCE_DIR}/tests/mock_archive.py stop ${TIAGO_MOCK_DIR})
set_tests_properties(tiago_mock_stop PROPERTIES FIXTURES_CLEANUP archive_http)
add_executable(test_http_client tests/test_http_client.f90)
target_link_libraries(test_http_client PRIVATE tiago_signals)
add_test(NAME tiago_http_client COMMAND test_http_client ${TIAGO_MOCK_DIR})
add_executable(test_online_source tests/test_online_source.f90)
target_link_libraries(test_online_source PRIVATE tiago_signals)
add_test(NAME tiago_online_source COMMAND test_online_source ${TIAGO_MOCK_DIR})
set_tests_properties(tiago_http_client tiago_online_source PROPERTIES
    FIXTURES_REQUIRED archive_http TIMEOUT 60 RESOURCE_LOCK archive_http)
set_tests_properties(tiago_http_client PROPERTIES FIXTURES_SETUP archive_client)
set_tests_properties(tiago_online_source PROPERTIES
    FIXTURES_REQUIRED "archive_http;archive_client" FIXTURES_SETUP archive_online)
add_executable(test_http_parallel tests/test_http_parallel.f90)
target_link_libraries(test_http_parallel PRIVATE tiago_signals)
add_test(NAME tiago_http_parallel COMMAND test_http_parallel ${TIAGO_MOCK_DIR})
set_tests_properties(tiago_http_parallel PROPERTIES FIXTURES_REQUIRED "archive_http;archive_online"
    FIXTURES_SETUP archive_parallel
    DEPENDS tiago_online_source TIMEOUT 30 RESOURCE_LOCK archive_http)
add_test(NAME tiago_archive_events COMMAND ${Python3_EXECUTABLE}
    ${TIAGO_SOURCE_DIR}/tests/check_archive_events.py ${TIAGO_MOCK_DIR})
set_tests_properties(tiago_archive_events PROPERTIES FIXTURES_REQUIRED "archive_http;archive_parallel"
    FIXTURES_SETUP archive_checked
    DEPENDS "tiago_http_client;tiago_online_source;tiago_http_parallel")

add_executable(tiago_online_example examples/online_signals.f90)
target_link_libraries(tiago_online_example PRIVATE tiago_inference)
add_test(NAME tiago_online_example COMMAND tiago_online_example
    ${TIAGO_MOCK_DIR}/adapter.json)
set_tests_properties(tiago_online_example PROPERTIES FIXTURES_REQUIRED "archive_http;archive_checked"
    DEPENDS tiago_archive_events TIMEOUT 30 RESOURCE_LOCK archive_http)

add_executable(benchmark_archive EXCLUDE_FROM_ALL
    benchmarks/signals/benchmark_archive.f90)
target_link_libraries(benchmark_archive PRIVATE tiago_signals)

add_executable(test_ece_transfer tests/test_ece_transfer.f90)
target_link_libraries(test_ece_transfer PRIVATE tiago_inference)
add_test(NAME tiago_ece_transfer COMMAND test_ece_transfer)
add_executable(benchmark_ece_transfer EXCLUDE_FROM_ALL benchmarks/ece/benchmark_transfer.f90)
target_link_libraries(benchmark_ece_transfer PRIVATE tiago_inference)
add_executable(benchmark_ece_receiver EXCLUDE_FROM_ALL benchmarks/ece/benchmark_receiver.f90)
target_link_libraries(benchmark_ece_receiver PRIVATE tiago_inference)

add_library(tiago_ece_geometry_fixture STATIC tests/ece_geometry_fixture.f90)
target_link_libraries(tiago_ece_geometry_fixture PUBLIC tiago_inference)
add_executable(test_ece_geometry tests/test_ece_geometry.f90)
target_link_libraries(test_ece_geometry PRIVATE tiago_ece_geometry_fixture)
add_test(NAME tiago_ece_geometry COMMAND $<TARGET_FILE:test_ece_geometry>)
add_test(NAME tiago_ece_geometry_vmecpp COMMAND $<TARGET_FILE:test_ece_geometry> vmecpp)
add_test(NAME tiago_ece_geometry_vmec COMMAND $<TARGET_FILE:test_ece_geometry>
    ${TIAGO_WOUT_SIMSOPT_SMALL})
set_tests_properties(tiago_ece_geometry_vmec PROPERTIES FIXTURES_REQUIRED wout_li383)
if(FORTSYM_BUILD)
    add_test(NAME tiago_ece_geometry_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ece_geometry> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/ece-geometry-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_ece_geometry_codegen.cmake)
endif()

add_executable(test_ece_geometry_limits tests/test_ece_geometry_limits.f90)
target_link_libraries(test_ece_geometry_limits PRIVATE tiago_ece_geometry_fixture)
add_test(NAME tiago_ece_geometry_limits COMMAND $<TARGET_FILE:test_ece_geometry_limits>)

add_executable(benchmark_ece_geometry EXCLUDE_FROM_ALL benchmarks/ece/benchmark_geometry.f90)
target_link_libraries(benchmark_ece_geometry PRIVATE tiago_ece_geometry_fixture)

add_executable(test_ece_opacity tests/test_ece_opacity.f90)
target_link_libraries(test_ece_opacity PRIVATE tiago_inference)
add_test(NAME tiago_ece_opacity COMMAND test_ece_opacity
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/ece_x1_flr.dat)
if(FORTSYM_BUILD)
    add_test(NAME tiago_ece_opacity_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ece_opacity> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/ece-opacity-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_ece_opacity_codegen.cmake)
endif()
add_executable(benchmark_ece_opacity EXCLUDE_FROM_ALL benchmarks/ece/benchmark_opacity.f90)
target_link_libraries(benchmark_ece_opacity PRIVATE tiago_inference)

if(FORTSYM_BUILD)
    add_test(NAME tiago_ece_recorded_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_ece_recorded> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/ece-recorded-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_ece_recorded_codegen.cmake)
endif()

add_library(tiago_ece_recorded_fixture STATIC tests/ece_recorded_fixture.f90)
target_link_libraries(tiago_ece_recorded_fixture PUBLIC tiago_ece_geometry_fixture)
add_executable(test_ece_recorded tests/test_ece_recorded.f90)
target_link_libraries(test_ece_recorded PRIVATE tiago_ece_recorded_fixture)
add_test(NAME tiago_ece_recorded COMMAND test_ece_recorded)
add_executable(test_ece_posterior tests/test_ece_posterior.f90)
target_link_libraries(test_ece_posterior PRIVATE tiago_ece_recorded_fixture)
add_test(NAME tiago_ece_posterior COMMAND test_ece_posterior
    ${TIAGO_BINARY_DIR}/tests/ece-posterior.csv)
add_executable(benchmark_ece_recorded EXCLUDE_FROM_ALL benchmarks/ece/benchmark_recorded.f90)
target_link_libraries(benchmark_ece_recorded PRIVATE tiago_ece_recorded_fixture)

add_executable(test_robust_signals tests/test_robust_signals.f90)
target_link_libraries(test_robust_signals PRIVATE tiago_inference)
add_test(NAME tiago_robust_signals COMMAND test_robust_signals)
if(TARGET gen_robust_signals)
    add_test(NAME tiago_robust_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_robust_signals> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/robust-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_robust_codegen.cmake)
endif()

if(TIAGO_ENABLE_RECONSTRUCTION)
    add_library(tiago_joint_recorded_fixture STATIC
        tests/joint_recorded_fixture.f90 tests/joint_signal_io.f90
        tests/joint_validation.f90)
    target_link_libraries(tiago_joint_recorded_fixture PUBLIC tiago_reconstruction tiago_signals)
    add_executable(test_joint_recorded tests/test_joint_recorded.f90)
    target_link_libraries(test_joint_recorded PRIVATE tiago_joint_recorded_fixture)
    add_test(NAME tiago_joint_recorded COMMAND ${TIAGO_SOURCE_DIR}/tests/joint_recorded.sh
        $<TARGET_FILE:test_joint_recorded> ${TIAGO_INPUT_LI383}
        ${TIAGO_SOURCE_DIR}/benchmarks/reconstruction/sensors.awk
        ${TIAGO_BINARY_DIR}/tests/joint-recorded)
    set_tests_properties(tiago_joint_recorded PROPERTIES FIXTURES_REQUIRED input_li383 TIMEOUT 3600)
endif()

if(TIAGO_IS_TOP_LEVEL)
    FetchContent_GetProperties(fortgen_runtime)
    add_test(NAME tiago_subproject COMMAND ${CMAKE_COMMAND}
        -DSOURCE=${TIAGO_SOURCE_DIR} -DOUTPUT=${TIAGO_BINARY_DIR}/tests/subproject
        -DGENERATOR=${CMAKE_GENERATOR}
        -DC_COMPILER=${CMAKE_C_COMPILER} -DFORTRAN_COMPILER=${CMAKE_Fortran_COMPILER}
        -DBUILD_TYPE=${CMAKE_BUILD_TYPE} -DCONFIG=$<CONFIG>
        -DFORTNUM_SOURCE=${fortnum_SOURCE_DIR} -DLIBNEO_SOURCE=${libneo_SOURCE_DIR}
        -DFORTIO_SOURCE=${fortio_SOURCE_DIR} -DFORTGEN_SOURCE=${fortgen_runtime_SOURCE_DIR}
        -DJSONFORTRAN_SOURCE=${jsonfortran_SOURCE_DIR}
        -P ${TIAGO_SOURCE_DIR}/tests/check_subproject.cmake)
    set_tests_properties(tiago_subproject PROPERTIES TIMEOUT 600)
endif()

# Energetic-particle event/surface response is independent of profile diagnostics.
add_library(tiago_wall STATIC
    src/generated/wall_energy.f90
    src/generated/wall_flux.f90
    src/generated/wall_thermal.f90
    src/generated/wall_radiance.f90
    src/diagnostics/wall_synthetic.f90)
set_target_properties(tiago_wall PROPERTIES
    Fortran_STANDARD 2018 Fortran_STANDARD_REQUIRED YES
    Fortran_MODULE_DIRECTORY "${TIAGO_MODULE_DIR}")
target_include_directories(tiago_wall PUBLIC
    "$<BUILD_INTERFACE:${TIAGO_MODULE_DIR}>")
target_link_libraries(tiago_wall PUBLIC tiago_signals)
if(TIAGO_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
    target_link_libraries(tiago_wall PUBLIC OpenMP::OpenMP_Fortran)
endif()
add_executable(test_wall_synthetic tests/test_wall_synthetic.f90)
target_link_libraries(test_wall_synthetic PRIVATE tiago_wall)
add_test(NAME tiago_wall_synthetic COMMAND test_wall_synthetic)
set_tests_properties(tiago_wall_synthetic PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
add_executable(benchmark_wall EXCLUDE_FROM_ALL benchmarks/benchmark_wall.f90)
target_link_libraries(benchmark_wall PRIVATE tiago_wall)
if(FORTSYM_BUILD)
    add_test(NAME tiago_wall_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_wall> -DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}
        -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/tests/wall-codegen
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_wall_codegen.cmake)
endif()
target_sources(tiago_wall PRIVATE src/generated/wall_legacy_energy.f90
    src/diagnostics/simple_wall_source.f90)
target_link_libraries(tiago_wall PRIVATE ${TIAGO_NETCDF_C_LIBRARY})
add_executable(test_simple_wall_source tests/test_simple_wall_source.f90)
target_link_libraries(test_simple_wall_source PRIVATE tiago_wall)
add_test(NAME tiago_simple_wall_fixture COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/simple_wall_fixture.py ${CMAKE_CURRENT_BINARY_DIR}/tests/simple-wall)
set_tests_properties(tiago_simple_wall_fixture PROPERTIES FIXTURES_SETUP simple_wall_file)
add_test(NAME tiago_simple_wall_source COMMAND test_simple_wall_source
    ${CMAKE_CURRENT_BINARY_DIR}/tests/simple-wall)
set_tests_properties(tiago_simple_wall_source PROPERTIES FIXTURES_REQUIRED simple_wall_file)

add_executable(test_magnetic_pickup tests/test_magnetic_pickup.f90)
target_link_libraries(test_magnetic_pickup PRIVATE tiago_inference)
add_test(NAME tiago_magnetic_pickup COMMAND test_magnetic_pickup)

add_executable(benchmark_magnetic_pickup EXCLUDE_FROM_ALL
    benchmarks/benchmark_magnetic_pickup.f90)
target_link_libraries(benchmark_magnetic_pickup PRIVATE tiago_inference)

if(FORTSYM_BUILD)
    add_test(NAME tiago_magnetic_pickup_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_magnetic_pickup> -DSOURCE=${TIAGO_SOURCE_DIR}
        -DOUTPUT=${TIAGO_BINARY_DIR}/tests/magnetic-pickup-codegen
        -P ${TIAGO_SOURCE_DIR}/tests/check_magnetic_pickup_codegen.cmake)
endif()

add_executable(test_harmonic_lockin tests/test_harmonic_lockin.f90)
target_link_libraries(test_harmonic_lockin PRIVATE tiago_inference)
add_test(NAME tiago_harmonic_lockin COMMAND test_harmonic_lockin)
if(FORTSYM_BUILD)
    add_test(NAME tiago_harmonic_codegen COMMAND ${CMAKE_COMMAND}
        -DGENERATOR=$<TARGET_FILE:gen_harmonic_lockin> -DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}
        -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/tests/harmonic-codegen
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_harmonic_codegen.cmake)
endif()
