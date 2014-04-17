## -------------------------------------------------------------------
## Copyright 2013 qdove.
##
## Author: Toby D. Young
## -------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

# qdove
set (QDOVE_VERSION 0.3.0)
set (QDOVE_BASE_NAME qdove)
set (CMAKE_BUILD_TYPE Debug)

include (CheckFunctionExists)
include (CheckSymbolExists)
include (CheckCXXCompilerFlag)
include (CheckCXXSourceCompiles)
include (CheckCXXSourceRuns)
include (CheckIncludeFile)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)

include (FindPackageHandleStandardArgs)

## doxygen and friends
find_package (Doxygen)
if (DOXYGEN_FOUND)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doc/qdove.dox.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/qdove.dox @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/qdove.dox
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with doxygen" VERBATIM
    )
endif (DOXYGEN_FOUND)

## deal.II 8.0 or newer
find_package (deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )


# Are all deal.II dependencies fullfilled?
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_SLEPC)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_SLEPC = ON
One or all of these are OFF in your installation but are required for this tutorial step."
    )
ENDIF()


DEAL_II_INITIALIZE_CACHED_VARIABLES ()
include_directories (${DEAL_II_INCLUDE_DIRS})

set (CMAKE_INCLUDE_CURRENT_DIR true)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set (QDOVE_WITH_LIBRARIES ${QDOVE_WITH_LIBRARIES} ${DEAL_II_LIBRARIES})

# start project
project (QDOVE)

enable_language (C)
enable_language (CXX)

set (LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
add_subdirectory (source)

# Appendix


