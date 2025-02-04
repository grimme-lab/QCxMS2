# This file is part of QCxMS2.

#########################################################################################
#########################################################################################
# Add modules to the CMake build 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

# specify module installation directory
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Options for enabling or disabling features
option(WITH_OpenMP  "Enable OpenMP support"       TRUE)
option(WITH_TOMLF   "Enable support for toml-f"   TRUE)
option(WITH_MCTCLIB   "Enable support for mctc-lib"   TRUE)

option(STATICBUILD  "Attempt to link everything statically" FALSE) # doesn't work yet

#########################################################################################
#########################################################################################

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

#########################################################################################
#########################################################################################

# Compiler settings for GNU and Intel Fortran compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check -ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-g -O2 -r8 -align array64byte -traceback")
  set(bounds "-check all -fpe0")
else()
  message(FATAL_ERROR "Please use an Intel or GNU compiler!")
endif()

# Apply the compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)

#########################################################################################
#########################################################################################

# Populate metadata.fh
set(version ${PROJECT_VERSION}${SOVERSION})
execute_process(COMMAND git show -s --format=%h
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(git_return)
  set(commit "unknown-commit")
endif()
string(TIMESTAMP date "%a, %d %B %H:%M:%S, %m/%d/%Y")
set(author $ENV{USERNAME})
execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
set(origin ${HOSTNAME})
set(fcid ${CMAKE_Fortran_COMPILER_ID})
set(fcver ${CMAKE_Fortran_COMPILER_VERSION})
set(ccid ${CMAKE_C_COMPILER_ID})
set(ccver ${CMAKE_C_COMPILER_VERSION})
set(bsystem "cmake ${CMAKE_VERSION}")
set(tomlfvar "${WITH_TOMLF}")
set(mctclibvar "${WITH_MCTCLIB}")

configure_file(
  "${PROJECT_SOURCE_DIR}/assets/templates/metadata.f90"
  "${PROJECT_BINARY_DIR}/qcxms2_version.fh"
  @ONLY
)
## just to be safe, create it also in include/
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/templates/metadata.f90"
  "${PROJECT_BINARY_DIR}/include/qcxms2_version.fh"
  @ONLY
)

#########################################################################################
#########################################################################################

