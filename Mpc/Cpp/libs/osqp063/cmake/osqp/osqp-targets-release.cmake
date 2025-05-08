#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "osqp::osqpstatic" for configuration "Release"
set_property(TARGET osqp::osqpstatic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(osqp::osqpstatic PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/osqp.lib"
  )

list(APPEND _cmake_import_check_targets osqp::osqpstatic )
list(APPEND _cmake_import_check_files_for_osqp::osqpstatic "${_IMPORT_PREFIX}/lib/osqp.lib" )

# Import target "osqp::osqp" for configuration "Release"
set_property(TARGET osqp::osqp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(osqp::osqp PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/osqp.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/osqp.dll"
  )

list(APPEND _cmake_import_check_targets osqp::osqp )
list(APPEND _cmake_import_check_files_for_osqp::osqp "${_IMPORT_PREFIX}/lib/osqp.lib" "${_IMPORT_PREFIX}/bin/osqp.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
