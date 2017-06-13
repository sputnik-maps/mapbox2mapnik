#
# Finds Mapbox
#
#  Mapbox_INCLUDE_DIR - Where to find Mapbox headers
#  Mapbox_FOUND       - True if Mapbox found
#

FIND_PATH(Mapbox_INCLUDE_DIR mapbox/variant.hpp PATH_SUFFIXES mapnik)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mapbox
    FOUND_VAR Mapbox_FOUND
    REQUIRED_VARS
        Mapbox_INCLUDE_DIR
)

MARK_AS_ADVANCED(Mapbox_INCLUDE_DIR)
