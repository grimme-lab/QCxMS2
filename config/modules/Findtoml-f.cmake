set(_lib "toml-f")
set(_pkg "TOML-F")
set(_url "https://github.com/toml-f/toml-f")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/qcxms2-utils.cmake")

qcxms2_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

if(TARGET "toml-f::toml-f")
  set (found TRUE)
else()
  set (found FALSE)
endif()
message(STATUS "Found toml-f: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
