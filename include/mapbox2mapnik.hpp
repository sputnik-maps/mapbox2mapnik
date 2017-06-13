#pragma once

#include <string>

#include <mapnik/map.hpp>

#include <mapbox2mapnik/config.hpp>


namespace sputnik {

void MAPBOX2MAPNIK_EXPORT load_mapbox_map(mapnik::Map& map, std::string const& filename,
                                          bool strict = false, const std::string& base_path="");

void MAPBOX2MAPNIK_EXPORT load_mapbox_map_string(mapnik::Map& map, std::string const& str,
                                                 bool strict = false, const std::string& base_path="");

} // ns sputnik
