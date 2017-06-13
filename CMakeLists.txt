cmake_minimum_required(VERSION 3.0.2)

project(mapbox2mapnik)

OPTION(ENABLE_SIMPLE_PLACEMENTS "Enable support for simple placements" OFF)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Boost MODULE COMPONENTS filesystem REQUIRED)
find_package(Glog MODULE REQUIRED)
find_package(Jsoncpp MODULE REQUIRED)
find_package(ICU REQUIRED)
find_package(Mapbox REQUIRED)
find_package(Mapnik 3.0 MODULE REQUIRED)

include(GenerateExportHeader)
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif (CMAKE_COMPILER_IS_GNUCXX)
if (ENABLE_SIMPLE_PLACEMENTS)
    add_definitions("-DMAPBOX2MAPNIK_ENABLE_SIMPLE_PLACEMENTS")
endif(ENABLE_SIMPLE_PLACEMENTS)
add_definitions("${Mapnik_DEFINITIONS}")
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(mapbox2mapnik_src src/mapbox2mapnik.cpp include/mapbox2mapnik.hpp)
add_library(mapbox2mapnik SHARED ${mapbox2mapnik_src})
set_property(TARGET mapbox2mapnik PROPERTY CXX_STANDARD 11)
set_property(TARGET mapbox2mapnik PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(mapbox2mapnik PUBLIC
    ${CMAKE_SOURCE_DIR}/include/
    ${PROJECT_BINARY_DIR}
    ${Boost_INCLUDE_DIRS}
    ${Glog_INCLUDE_DIR}
    ${Jsoncpp_INCLUDE_DIR}
    ${ICU_INCLUDE_DIR}
    ${Mapbox_INCLUDE_DIR}
    ${Mapnik_INCLUDE_DIRS})
target_link_libraries(mapbox2mapnik
    ${Mapnik_LIBRARIES}
    ${Glog_LIBRARY}
    ${ICU_LIBRARY}
    ${Boost_LIBRARIES}
    ${Jsoncpp_LIBRARIES})
generate_export_header(mapbox2mapnik EXPORT_FILE_NAME mapbox2mapnik/config.hpp)

install (TARGETS mapbox2mapnik DESTINATION lib)
install (FILES
    include/mapbox2mapnik.hpp
    ${PROJECT_BINARY_DIR}/mapbox2mapnik/config.hpp DESTINATION include/mapbox2mapnik)
