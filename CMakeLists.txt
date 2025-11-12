cmake_minimum_required(VERSION 3.24)
project(tiago LANGUAGES Fortran)

include(FetchContent)

set(TIAGO_LIBNEO_GIT "git@github.com:itpplasma/libneo.git" CACHE STRING "libneo upstream to fetch")
set(TIAGO_LIBNEO_TAG "main" CACHE STRING "libneo branch/tag used by Tiago")
set(LIBNEO_ENABLE_TESTS OFF CACHE BOOL "Disable libneo's own tests when vendored" FORCE)
set(ENABLE_OPENMP OFF CACHE BOOL "Disable libneo OpenMP flags unless explicitly requested" FORCE)

FetchContent_Declare(libneo
    GIT_REPOSITORY ${TIAGO_LIBNEO_GIT}
    GIT_TAG        ${TIAGO_LIBNEO_TAG}
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  .)
FetchContent_MakeAvailable(libneo)

add_library(tiago_core INTERFACE)
target_link_libraries(tiago_core INTERFACE LIBNEO::neo)

enable_testing()
add_custom_target(tiago_tests ALL
    COMMENT "Tests will be added in Stage 2 once the vacuum solver lands")
