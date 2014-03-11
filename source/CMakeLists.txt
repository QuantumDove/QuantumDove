
## Main trunk of gubbins includes
include_directories (
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}
  )

## Subdirectoryies in the source tree
add_subdirectory (base)
add_subdirectory (generic_linear_algebra)
add_subdirectory (materials)
add_subdirectory (models)
add_subdirectory (psuedopotentials)

## Wrap objects into a single library
add_library (${GUBBINS_BASE_NAME} SHARED ## NB: this "SHARED" should not really be here!
  $<TARGET_OBJECTS:base>
  $<TARGET_OBJECTS:generic_linear_algebra>
  $<TARGET_OBJECTS:materials>
  $<TARGET_OBJECTS:models>
  $<TARGET_OBJECTS:psuedopotentials>
)

## Set library properties
set_target_properties (${GUBBINS_BASE_NAME} PROPERTIES SOVERSION ${GUBBINS_VERSION})

## Link (with external) libraries
# target_link_libraries (${GUBBINS_BASE_NAME} ${GUBBINS_EXTERNAL_LIBRARIES})
