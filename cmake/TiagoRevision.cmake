# Writes OUTPUT with the source revision; rewritten only when it changes so
# that an unchanged revision triggers no rebuild.
execute_process(COMMAND git describe --always --dirty --abbrev=40
    WORKING_DIRECTORY "${SOURCE_DIR}"
    OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE status ERROR_QUIET)
if(NOT status EQUAL 0 OR revision STREQUAL "")
    set(revision "unknown")
endif()
set(content "module tiago_revision
    implicit none
    character(len=*), parameter :: tiago_source_revision = \"${revision}\"
end module tiago_revision
")
if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" old)
else()
    set(old "")
endif()
if(NOT old STREQUAL content)
    file(WRITE "${OUTPUT}" "${content}")
endif()
