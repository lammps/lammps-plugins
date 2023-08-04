##########################################
# CMake build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

project(lammpsplugins VERSION 1.0 LANGUAGES CXX)

if(NOT LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Must set LAMMPS_SOURCE_DIR variable")
endif()
set(CMAKE_MODULE_PATH "${LAMMPS_SOURCE_DIR}/../cmake/Modules")
include(CheckIncludeFileCXX)
include(LAMMPSInterfacePlugin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  add_custom_target(package DEPENDS package-USER-AEAM package-USER-BFIELD package-USER-REBOMOS)
endif()
add_subdirectory(USER-AEAM)
add_subdirectory(USER-BFIELD)
add_subdirectory(USER-REBOMOS)
