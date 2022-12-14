###################################################################################
#
# Vessel tree construction
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

cmake_minimum_required(VERSION 3.16)

option(enable_OPENMP "Enable OpenMP" ON)
option(enable_GLOG "Enable GLOG" OFF)

#-----compiler options
if(NOT DEFINED compiler)
  set(base_compile_options "-Wall -Wextra -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -fopenmp")
   endif()
endif()

if(compiler STREQUAL intel)
  set(CMAKE_C_COMPILER "icx")
  set(CMAKE_CXX_COMPILER "icpx")
  set(base_compile_options "-Wall -Wextra -O3")
  #if(enable_OPENMP)
  set(base_compile_options ${base_compile_options} " -qopenmp")
  set(base_compile_options ${base_compile_options} " -pg")
  #endif()
elseif(DEFINED compiler)
  message(FATAL ERROR "invalid compiler")
endif()
string(REPLACE " " " " compile_options ${base_compile_options})
#-----------------------

project(fem_homozinied_diffusion C CXX)
set(CMAKE_CXX_STANDARD 11)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

#change default install prefix
set(INSTALL_DIR "/mnt/d/diffusion")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "..." FORCE)
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

#######
# Display options
#######
message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

# TextParser
set(TP_DIR ${TP_DIR} CACHE PATH "...")
find_package(TP REQUIRED)
include_directories(${TP_INCLUDE_DIR})

# HDF5
find_package (HDF5 REQUIRED COMPONENTS C CXX)
IF(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  MESSAGE("find HDF5")
  MESSAGE(${HDF5_INCLUDE_DIRS})
  MESSAGE(${HDF5_LIBRARIES})
  MESSAGE(${HDF5_CXX_LIBRARIES})
ELSE()
  MESSAGE(FATAL_ERROR "not find HDF5")
ENDIF(HDF5_FOUND)

# MKL
set(MKL_ARCH "intel64")
set(MKL_LINK "static")
set(MKL_INTERFACE "lp64")
set(MKL_THREADING "intel_thread")
set(MKL_MPI "intelmpi")
find_package(MKL REQUIRED C CXX)
include_directories(${MKL_INCLUDE_DIR})
MESSAGE(${MKL_LIBRARIES})

#######
# Sub Directories
#######

# add_subdirectory(include)
add_subdirectory(src)