include_guard(GLOBAL)

# External equilibrium providers reuse these targets without the VMEC stack.
if(NOT TARGET fortnum)
    message(FATAL_ERROR "TiagoInferenceCore requires the fortnum target")
endif()
find_package(LAPACK REQUIRED)
set(_tiago_core_source "${CMAKE_CURRENT_LIST_DIR}/../src")

add_library(tiago_inference_core STATIC
    "${_tiago_core_source}/inference/local_posterior.f90"
    "${_tiago_core_source}/inference/linear_posterior.f90"
    "${_tiago_core_source}/inference/signal_correlation.f90"
    "${_tiago_core_source}/inference/inference_engines.f90")
target_link_libraries(tiago_inference_core PUBLIC fortnum ${LAPACK_LIBRARIES})

add_library(tiago_point_diagnostics STATIC
    "${_tiago_core_source}/diagnostics/diagnostic_types.f90"
    "${_tiago_core_source}/generated/mse_rows.f90"
    "${_tiago_core_source}/generated/mse_local.f90"
    "${_tiago_core_source}/diagnostics/mse.f90")

foreach(target IN ITEMS tiago_inference_core tiago_point_diagnostics)
    set(_tiago_core_modules "${CMAKE_CURRENT_BINARY_DIR}/${target}-modules")
    file(MAKE_DIRECTORY "${_tiago_core_modules}")
    set_target_properties(${target} PROPERTIES
        Fortran_STANDARD 2018
        Fortran_STANDARD_REQUIRED YES
        Fortran_MODULE_DIRECTORY "${_tiago_core_modules}")
    target_include_directories(${target} PUBLIC
        "$<BUILD_INTERFACE:${_tiago_core_modules}>")
endforeach()
