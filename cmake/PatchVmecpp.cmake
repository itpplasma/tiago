# Applied to the fetched VMEC++ sources (FetchContent PATCH_COMMAND):
#  - drop the pybind11 extension module, so that no Python is configured or
#    built; Tiago uses the C++ core through src/equilibrium/vmecpp_adapter.cc;
#  - optionally replace the HDF5 and netCDF-C source archives by mirrors
#    (-DHDF5_URL=..., -DNETCDF_URL=...; the upstream hashes are then dropped).
# Idempotent: a second application finds nothing left to change.
file(READ "${SOURCE}/CMakeLists.txt" text)

set(begin_marker "# Now add the pybind11 module for VMEC++.")
set(end_marker "install(TARGETS _vmecpp LIBRARY DESTINATION vmecpp/cpp/.)")
string(FIND "${text}" "${begin_marker}" begin)
string(FIND "${text}" "${end_marker}" end)
if(begin GREATER -1 AND end GREATER begin)
    string(LENGTH "${end_marker}" end_length)
    math(EXPR end "${end} + ${end_length}")
    string(SUBSTRING "${text}" 0 ${begin} head)
    string(SUBSTRING "${text}" ${end} -1 tail)
    set(text "${head}# (pybind11 module removed by Tiago's PatchVmecpp.cmake)${tail}")
endif()

if(HDF5_URL)
    string(REGEX REPLACE "URL \"https://github.com/HDFGroup/hdf5/archive/refs/tags/[^\"]*\""
        "URL \"${HDF5_URL}\"" text "${text}")
    string(REGEX REPLACE "\n[ \t]*URL_HASH SHA256=df5ee33c74d5efb59738075ef96f4201588e1f1eeb233f047ac7fd1072dee1f6"
        "" text "${text}")
endif()
if(NETCDF_URL)
    string(REGEX REPLACE "URL \"https://github.com/Unidata/netcdf-c/archive/refs/tags/[^\"]*\""
        "URL \"${NETCDF_URL}\"" text "${text}")
    string(REGEX REPLACE "\n[ \t]*URL_HASH SHA256=990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff"
        "" text "${text}")
endif()

file(WRITE "${SOURCE}/CMakeLists.txt" "${text}")
