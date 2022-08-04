##########################################
# CMake build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

project(aeamplugin VERSION 1.0 LANGUAGES CXX)

if(NOT LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Must set LAMMPS_SOURCE_DIR variable")
endif()
set(CMAKE_MODULE_PATH "${LAMMPS_SOURCE_DIR}/../cmake/Modules")
if (NOT TARGET lammps)
  include(CheckIncludeFileCXX)
  include(LAMMPSInterfacePlugin)
endif()

##########################
# building the plugins

add_library(aeamplugin MODULE aeamplugin.cpp ${CMAKE_CURRENT_SOURCE_DIR}/pair_aeam.cpp)
target_link_libraries(aeamplugin PRIVATE lammps)
set_target_properties(aeamplugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(aeamplugin PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(aeamplugin PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(aeamplugin  PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()

  get_lammps_version(${LAMMPS_SOURCE_DIR}/version.h LAMMPS_VERSION)
  find_program(MAKENSIS_PATH makensis)
  if(MAKENSIS_PATH)
      execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/../lammps.ico
      ${CMAKE_CURRENT_SOURCE_DIR}/../lammps-text-logo-wide.bmp ${CMAKE_CURRENT_SOURCE_DIR}/aeamplugin.nsis ${CMAKE_BINARY_DIR})
    if(TARGET package)
      set(mytarget package-USER-AEAM)
    else()
      set(mytarget package)
    endif()
    if(BUILD_MPI)
      add_custom_target(${mytarget} ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MPI aeamplugin.nsis
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DEPENDS aeamplugin ${CMAKE_BINARY_DIR}/lammps.ico ${CMAKE_BINARY_DIR}/lammps-text-logo-wide.bmp ${CMAKE_BINARY_DIR}/aeamplugin.nsis
        BYPRODUCTS ${CMAKE_BINARY_DIR}/LAMMPS-USER-AEAM-plugin-${LAMMPS_VERSION}-MPI.exe)
    else()
      add_custom_target(${mytarget} ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION} aeamplugin.nsis
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DEPENDS aeamplugin ${CMAKE_BINARY_DIR}/lammps.ico ${CMAKE_BINARY_DIR}/lammps-text-logo-wide.bmp ${CMAKE_BINARY_DIR}/aeamplugin.nsis
        BYPRODUCTS ${CMAKE_BINARY_DIR}/LAMMPS-USER-AEAM-plugin-${LAMMPS_VERSION}.exe)
    endif()
  endif()
else()
  set_target_properties(aeamplugin PROPERTIES LINK_FLAGS "-rdynamic")
endif()
