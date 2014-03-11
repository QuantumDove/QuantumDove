## -------------------------------------------------------------------
## Copyright 2013 gubbins.
##
## Author: Toby D. Young
## -------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

# gubbins
set (GUBBINS_VERSION 0.3.0)
set (GUBBINS_BASE_NAME qdove)
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
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/doc/gubbins.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/gubbins.dox @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doc/gubbins.dox
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
    COMMENT "Generating API documentation with doxygen" VERBATIM
    )
endif (DOXYGEN_FOUND)

## deal.II 8.0 or newer
find_package (deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES ()
include_directories (${DEAL_II_INCLUDE_DIRS})

set (CMAKE_INCLUDE_CURRENT_DIR true)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
set (GUBBINS_WITH_LIBRARIES ${GUBBINS_WITH_LIBRARIES} ${DEAL_II_LIBRARIES})

# start project
project (GUBBINS)

enable_language (C)
enable_language (CXX)

set (LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
add_subdirectory (source)

# Appendix

#if (CXX_COMPILER_FLAG_RDYNAMIC)
#  add_definitions ("-rdynamic")
#endif ()

#foreach (library ${GUBBINS_WITH_LIBRARIES})
#  target_link_libraries (${GUBBINS_BASE_NAME} ${library}) 
#  message (STATUS "Added library ${library}")
#endforeach ()

