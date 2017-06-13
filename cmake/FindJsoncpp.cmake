#
# Finds Jsoncpp
#
#  Jsoncpp_INCLUDE_DIRS   - Where to find Jsoncpp headers
#  Jsoncpp_LIBRARIES      - The path of the Jsoncpp library
#  Jsoncpp_FOUND          - True if Jsoncpp found
#


find_path(Jsoncpp_INCLUDE_DIR json/json.h PATH_SUFFIXES jsoncpp)

find_library(Jsoncpp_LIBRARY jsoncpp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Jsoncpp DEFAULT_MSG Jsoncpp_INCLUDE_DIR Jsoncpp_LIBRARY)

if (Jsoncpp_FOUND)
    set(Jsoncpp_INCLUDE_DIRS ${Jsoncpp_INCLUDE_DIR})
    set(Jsoncpp_LIBRARIES ${Jsoncpp_LIBRARY})
endif (Jsoncpp_FOUND)

mark_as_advanced (Jsoncpp_INCLUDE_DIR Jsoncpp_LIBRARY)
