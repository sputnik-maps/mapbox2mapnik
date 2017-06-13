#
#  Finds Mapnik
#
#  Mapnik_FOUND        - True if Mapnik found
#  Mapnik_INCLUDE_DIRS - Where to find Mapnik headers
#  Mapnik_LIBRARIES    - The path of the Mapnik library
#  Mapnik_DEFINITIONS  - Pre-processor defines for Mapnik build
#  Mapnik_VERSION      - The version of Mapnik library
#  Mapnik_PLUGINDIR    - Where to find Mapnik input plugins
#


FIND_PATH(Mapnik_INCLUDE_DIR mapnik/config.hpp)

FIND_LIBRARY(Mapnik_LIBRARY mapnik)

FIND_PROGRAM(Mapnik_CONFIG mapnik-config)

if (Mapnik_CONFIG)
    execute_process(COMMAND mapnik-config --version OUTPUT_VARIABLE Mapnik_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND mapnik-config --input-plugins OUTPUT_VARIABLE Mapnik_PLUGINDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND mapnik-config --defines OUTPUT_VARIABLE Mapnik_DEFINITIONS
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif (Mapnik_CONFIG)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mapnik
    FOUND_VAR Mapnik_FOUND
    REQUIRED_VARS
        Mapnik_INCLUDE_DIR
        Mapnik_LIBRARY
    VERSION_VAR Mapnik_VERSION)

if(Mapnik_FOUND)
  set(Mapnik_LIBRARIES ${Mapnik_LIBRARY})
  set(Mapnik_INCLUDE_DIRS ${Mapnik_INCLUDE_DIR})
endif()

MARK_AS_ADVANCED(Mapnik_LIBRARY Mapnik_INCLUDE_DIR)
