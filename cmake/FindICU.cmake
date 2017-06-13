#
# Finds ICU
#
#  ICU_INCLUDE_DIR - Where to find ICU headers
#  ICU_LIBRARY     - The path of the ICU library
#  ICU_FOUND       - True if ICU found
#

FIND_PATH(ICU_INCLUDE_DIR unicode/uversion.h)

FIND_LIBRARY(ICU_LIBRARY icuuc)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ICU DEFAULT_MSG ICU_LIBRARY ICU_INCLUDE_DIR)

MARK_AS_ADVANCED(ICU_LIBRARY ICU_INCLUDE_DIR)
