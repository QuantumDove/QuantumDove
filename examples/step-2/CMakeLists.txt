cmake_minimum_required (VERSION 2.8.8)
include (FindPackageHandleStandardArgs)

set (TARGET "step-2")
set (TARGET_SRC
  step-2.cc
)

find_package (deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()
project (${TARGET})

# Find qdove libraries		
find_library (QDOVE_LIBRARIES
  NAMES qdove 
  PATHS "${PROJECT_SOURCE_DIR}/../../lib"
  )
find_package_handle_standard_args ("qdove libraries" REQUIRED_VARS QDOVE_LIBRARIES)

include_directories (${PROJECT_SOURCE_DIR}/../../include ${DEAL_II_INCLUDE_DIRS})

add_executable (${TARGET} ${TARGET_SRC})
target_link_libraries (${TARGET} ${DEAL_II_LIBRARIES} ${QDOVE_LIBRARIES})




