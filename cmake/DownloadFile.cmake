if(NOT DEFINED URL)
    message(FATAL_ERROR "URL is required")
endif()

if(NOT DEFINED DEST)
    message(FATAL_ERROR "DEST is required")
endif()

get_filename_component(_dest_dir "${DEST}" DIRECTORY)
if(_dest_dir STREQUAL "")
    set(_dest_dir ".")
endif()
file(MAKE_DIRECTORY "${_dest_dir}")

if(EXISTS "${DEST}")
    file(SIZE "${DEST}" _existing_size)
    if(_existing_size GREATER 0)
        message(STATUS "Using cached download: ${DEST}")
        return()
    else()
        message(WARNING "Cached file has zero size, redownloading: ${DEST}")
        file(REMOVE "${DEST}")
    endif()
endif()

message(STATUS "Downloading ${URL}")
file(DOWNLOAD
    "${URL}"
    "${DEST}"
    TIMEOUT 120
    STATUS _download_status
    SHOW_PROGRESS)

list(GET _download_status 0 _status_code)
list(GET _download_status 1 _status_msg)
if(NOT _status_code EQUAL 0)
    file(REMOVE "${DEST}")
    message(FATAL_ERROR "Download failed (${_status_code}): ${_status_msg}")
endif()

file(SIZE "${DEST}" _final_size)
if(_final_size EQUAL 0)
    file(REMOVE "${DEST}")
    message(FATAL_ERROR "Downloaded file is empty: ${DEST}")
endif()

message(STATUS "Saved to ${DEST} (${_final_size} bytes)")
