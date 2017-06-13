#
# Finds Glog
#
#  Glog_INCLUDE_DIR - Where to find Glog headers
#  Glog_LIBRARY     - The path of the Glog library
#  Glog_FOUND       - True if Glog found
#

FIND_PATH(Glog_INCLUDE_DIR glog/logging.h)

FIND_LIBRARY(Glog_LIBRARY glog)


INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Glog DEFAULT_MSG Glog_LIBRARY Glog_INCLUDE_DIR)

MARK_AS_ADVANCED(Glog_LIBRARY Glog_INCLUDE_DIR)
