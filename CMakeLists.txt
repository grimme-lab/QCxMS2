# This file is part of QCxMS2
# SPDX-Identifier: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.17)
option(WITH_OBJECT "To build using object library" TRUE)
option(INSTALL_MODULES "Install Fortran module files to include directory." FALSE)
option(WITH_TOMLF "Enable TOML-F support" TRUE)   

# Buggy CMake versions
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0 AND CMAKE_VERSION VERSION_LESS 3.28.0)
   set(WITH_OBJECT FALSE)
endif()

# Setup the QCxMS2 Project
project(
  qcxms2 
  LANGUAGES "C" "Fortran"
  VERSION  1.0.0
  DESCRIPTION "Program package for the quantum mechanical calculation of EI mass spectra using automated reaction network exploration "
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include further configurations
add_subdirectory("config")


###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################
# Check a specific CMake targets and execute the
# corresponding Find scripts (located in config/modules/)

# LAPACK, BLAS, OpenMP (usually via shared dependencies)
if(STATICBUILD)
 set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
 set(CMAKE_EXE_LINKER_FLAGS "-static")
 set(CMAKE_LINK_SEARCH_START_STATIC TRUE)
 set(CMAKE_LINK_SEARCH_END_STATIC TRUE)
 set(BLA_STATIC ON)
endif()
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
  get_target_property(OpenMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran INTERFACE_LINK_LIBRARIES)
  message(STATUS "OpenMP::OpenMP_Fortran is linking the following libraries:")
  foreach(lib ${OpenMP_Fortran_LIBRARIES})
    message(STATUS "${lib}")
  endforeach() 
endif()

# Check if we are using OpenBLAS (need a precompiler definition if yes)
if(LAPACK_LIBRARIES)
  string(FIND "${LAPACK_LIBRARIES}" "openblas" _openblas_in_lapack)

  if(NOT _openblas_in_lapack EQUAL -1)
    message(STATUS "libopenblas was found as part of LAPACK")
    add_compile_definitions(WITH_OPENBLAS)
  endif()
endif()

# Fortran unit test interface (also used by other subprojects)
if(NOT TARGET "test-drive::test-drive" AND WITH_TESTS)
   find_package("test-drive" REQUIRED)
endif()


if(WITH_TOMLF)                                                                                                                                                  
   find_package("toml-f" REQUIRED)                                                                                                                               
   add_compile_definitions(WITH_TOMLF)                                                                                                                           
 endif()   

# TOML-F (needs to be imported before tblite)
if(NOT TARGET "toml-f::toml-f" AND WITH_TOMLF)
   find_package("toml-f" REQUIRED)
   add_compile_definitions(WITH_TOMLF)
endif()


#mctc-lib 
if(NOT TARGET "mctc-lib::mctc-lib" AND WITH_MCTCLIB)
   find_package("mctc-lib" REQUIRED)
   add_compile_definitions(WITH_MCTCLIB)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "$<$<VERSION_LESS:${CMAKE_VERSION},3.20>:BLAS::BLAS>"
  "toml-f::toml-f"
  "mctc-lib::mctc-lib"
)

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)
add_subdirectory("src")


if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
########################## OBJECT LIBRARY #####################################
###############################################################################
if(WITH_OBJECT AND NOT STATICBUILD)
   
   add_library(
      "${PROJECT_NAME}-object"
      OBJECT
      ${srcs}
   )

   # customize object library
   set_target_properties(
      "${PROJECT_NAME}-object"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      POSITION_INDEPENDENT_CODE ON
   )

   # link object library conditionally against required
   target_link_libraries(
      "${PROJECT_NAME}-object"
      PUBLIC
      $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
      $<$<BOOL:${WITH_MCTCLIB}>:mctc-lib::mctc-lib>
      $<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>
   )

   # include directories
   target_include_directories(
      "${PROJECT_NAME}-object"
      PUBLIC
      ${QCxMS2-config-dir}
      ${PROJECT_BINARY_DIR}
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )

endif()

###############################################################################
############################### Static Library ################################
###############################################################################
if(WITH_OBJECT AND NOT STATICBUILD)
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )
else()
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      ${srcs}
   )
endif()

# workaround
set(LINK_OpenMP FALSE)
if(WITH_OpenMP AND NOT STATICBUILD)
  set(LINK_OpenMP TRUE)
endif()

#somehow needed ...                                                                                                                                             

set(TOMLF_LIB_PATH "${CMAKE_SOURCE_DIR}/subprojects/toml-f")                                                                                              
link_directories(${TOMLF_LIB_PATH})    

target_link_libraries(
   "lib-${PROJECT_NAME}-static"
   PUBLIC
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   $<$<BOOL:${LINK_OpenMP}>:OpenMP::OpenMP_Fortran>
   $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
   $<$<BOOL:${WITH_MCTCLIB}>:mctc-lib::mctc-lib>
   $<$<BOOL:${STATICBUILD}>:-static>
)


set_target_properties(
   "lib-${PROJECT_NAME}-static"
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
   ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
   POSITION_INDEPENDENT_CODE ON
   OUTPUT_NAME "${PROJECT_NAME}"
 )


target_include_directories(
  "lib-${PROJECT_NAME}-static"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)


####################################################################################
############################# Shared Library #######################################
####################################################################################
if (WITH_OBJECT AND NOT STATICBUILD) 
   add_library(
      "lib-${PROJECT_NAME}-shared"
      SHARED
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )

   target_link_libraries(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      $<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>
      $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
      $<$<BOOL:${WITH_MCTCLIB}>:mctc-lib::mctc-lib>
   )

   set_target_properties(
      "lib-${PROJECT_NAME}-shared"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
      OUTPUT_NAME "${PROJECT_NAME}"
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_VERSION_MAJOR}"
   )

   target_include_directories(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )
endif()

###############################################################################
############################### Executables ###################################
###############################################################################
add_executable(
  ${PROJECT_NAME}-exe
  ${prog}
)

target_link_libraries(
   ${PROJECT_NAME}-exe
   PRIVATE
   "lib-${PROJECT_NAME}-static"
   $<$<BOOL:${STATICBUILD}>:-static>
)

set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)

target_include_directories(
   ${PROJECT_NAME}-exe 
   PRIVATE 
   ${PROJECT_SOURCE_DIR}/include
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)


###############################################################################
########################### Install Instructions ##############################
###############################################################################

# Build output artifacts
if (WITH_OBJECT AND NOT STATICBUILD)
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
#      "lib-${PROJECT_NAME}-shared"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
else()
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

