#include "mapbox2mapnik.hpp"

#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <streambuf>
#include <typeinfo>
#include <vector>

#include <glog/logging.h>

#include <boost/filesystem/operations.hpp>
#include <boost/optional.hpp>

#include <json/json.h>

#include <mapnik/box2d.hpp>
#include <mapnik/color.hpp>
#include <mapnik/datasource.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/expression_node.hpp>
#include <mapnik/feature_type_style.hpp>
#include <mapnik/font_engine_freetype.hpp>
#include <mapnik/layer.hpp>
#include <mapnik/map.hpp>
#include <mapnik/rule.hpp>
#include <mapnik/symbolizer.hpp>
#include <mapnik/symbolizer_enumerations.hpp>
#include <mapnik/symbolizer_keys.hpp>
#include <mapnik/transform_expression.hpp>

#include <mapnik/text/formatting/text.hpp>
#include <mapnik/text/placements/dummy.hpp>
#ifdef MAPBOX2MAPNIK_ENABLE_SIMPLE_PLACEMENTS
#include <mapnik/text/placements/simple.hpp>
#endif
#include <mapnik/text/text_layout.hpp>

#include <mapnik/util/fs.hpp>
#include <mapnik/util/variant.hpp>


namespace sputnik {

const uint kMaxZoom = 19;
const uint kEmsInPixel = 16;
const double kDefaultSpacing = 350.0;


static double zoom2scale(double zoom){
    return 559082264. / std::pow(2, zoom);
}

static double mapbox_zoom2scale(double zoom) {
    return zoom2scale(zoom + 1) + 1;
}

static constexpr double ems2pixels(double ems, double text_size) {
    return ems * text_size;
}

static mapnik::expr_node str2expr_node(const std::string& str) {
    static const mapnik::transcoder tr("utf-8");
    return tr.transcode(str.data(), str.length());
}

static inline boost::optional<Json::Value> InterpolateVal(const Json::Value& val1, const Json::Value& val2,
                                                          double i_factor) {
    // TODO: other types
    assert(i_factor > 0.0 && i_factor < 1.0);
    if (val1.isNumeric()) {
        if (!val2.isNumeric()) {
            LOG(ERROR) << "Stop values have different types " << val1 << ' ' << val2;
            return boost::none;
        }
        double v1 = val1.asDouble();
        double v2 = val2.asDouble();
        double iv = v1 + (v2 - v1) * i_factor;
        return Json::Value(iv);
    }
    return boost::none;
}

static inline double InterpolationFactor(double zoom_d, double pos, double i_base) {
    if (i_base == 1.0) {
        return pos / zoom_d;
    }
    return (std::pow(i_base, pos) - 1) / (std::pow(i_base, zoom_d) - 1);
}

void Interpolate(std::map<double, Json::Value>& zoom_props_map, double i_base,
                 double minzoom, double maxzoom) {
    for (auto zoom_props_iter = zoom_props_map.begin(), next_iter = std::next(zoom_props_iter);
         next_iter != zoom_props_map.end(); zoom_props_iter = next_iter++) {
        double s_zoom = zoom_props_iter->first;
        double e_zoom = next_iter->first;
        double zoom_d = e_zoom - s_zoom;
        for (double pos = 1.0; pos < zoom_d; pos += 1.0) {
            double i_factor = InterpolationFactor(zoom_d, pos, i_base);
            boost::optional<Json::Value> i_val = InterpolateVal(zoom_props_iter->second, next_iter->second,
                                                                i_factor);
            if (i_val) {
                zoom_props_map[s_zoom + pos] = std::move(*i_val);
            }
        }
    }
}

class PropertiesView {
public:
    explicit PropertiesView(const Json::Value& layer_properties) {
        minzoom_ = layer_properties.get("minzoom", 1).asDouble();
        maxzoom_ = layer_properties.get("maxzoom", 20).asDouble();
        const Json::Value& mapbox_paint = layer_properties["paint"];
        const Json::Value& mapbox_layout = layer_properties["layout"];
        const Json::Value& metadata = layer_properties["metadata"];
        layer_id_ = layer_properties["id"].asString();
        view_data_[minzoom_];
        extract_properties(mapbox_paint);
        extract_properties(mapbox_layout);
        extract_properties(metadata);

        if (props_data_.empty()) {
            valid_ = false;
            return;
        }

        // Populate view_data_
        for (auto& zoom_props : view_data_) {
            double zoom = zoom_props.first;
            auto& props = zoom_props.second;
            for (auto& prop_stops : props_data_) {
                auto prop_stops_itr = prop_stops.lower_bound(zoom);
                if (prop_stops_itr == prop_stops.end()) {
                    prop_stops_itr = prop_stops.begin();
                }
                props.push_back(&prop_stops_itr->second);
            }
        }
        iter_ = view_data_.end();
        valid_ = true;
    }

    inline bool next() {
        if (!valid_) {
            return false;
        }
        if (iter_ == view_data_.end()) {
            reset_cursor();
            return true;
        }
        if (std::next(iter_) == view_data_.end()) {
            return false;
        }
        ++iter_;
        return true;
    }

    inline void reset_cursor() {
        iter_ = view_data_.begin();
    }

    inline double minzoom() const {
        if (!valid_ || iter_ == view_data_.end()) {
            return 0.0;
        }
        return iter_->first;
    }

    inline double maxzoom() const {
        if (!valid_ || iter_ == view_data_.end()) {
            return 0.0;
        }
        auto next_iter = std::next(iter_);
        if (next_iter == view_data_.end()) {
            return maxzoom_;
        }
        return next_iter->first;
    }

    const Json::Value& property(const std::string& prop_name) const {
        if (!valid_ || iter_ == view_data_.end()) {
            return empty_val_;
        }
        auto names_count = prop_names_.size();
        for (uint i = 0; i < names_count; ++i) {
            if (prop_names_[i] == prop_name) {
                return *iter_->second[i];
            }
        }
        return empty_val_;
    }

    inline bool valid() const {
        return valid_;
    }

private:
    void extract_properties(const Json::Value& json_props) {
        for (auto prop_itr = json_props.begin(); prop_itr != json_props.end() ; ++prop_itr) {
            std::map<double, Json::Value> props_map;
            const Json::Value& prop_val = *prop_itr;
            if (prop_val.isObject() && prop_val.isMember("stops")) {
                const Json::Value& stops = prop_val["stops"];
                if (!stops.isArray() || stops.empty()) {
                    LOG(ERROR) << "Invalid stops!\n" << stops;
                    continue;
                }

                uint num_stops = stops.size();
                for (uint i = 0; i < num_stops; ++i) {
                    const Json::Value& stop_pair = stops[i];
                    if (stop_pair.size() != 2) {
                        LOG(ERROR) << "Invalid stop " << stop_pair << " in\n" << stops;
                        continue;
                    }
                    double zoom = stop_pair[0].asDouble();
                    if (zoom < minzoom_ || zoom >= maxzoom_) {
                        LOG(WARNING) << layer_id_ << " Stop " << zoom << " out of bounds (" <<
                                      minzoom_ << ", " << maxzoom_ << ')';
                    }
                    props_map[zoom] = stop_pair[1];
                }
                if (props_map.size() > 1) {
                    // Insert interpolated values
                    double i_base = prop_val.get("base", 1.0).asDouble();
                    Interpolate(props_map, i_base, minzoom_, maxzoom_);
                }
                // Insert zooms into view_data_
                for (auto props_map_iter = props_map.begin(); props_map_iter != props_map.end(); ++props_map_iter) {
                    double zoom = props_map_iter->first;
                    if (zoom >= minzoom_ && zoom < maxzoom_) {
                        view_data_[zoom];
                    }
                }
            } else {
                props_map[minzoom_] = prop_val;
            }
            if (props_map.empty()) {
                LOG(WARNING) << "No props found for layer " << layer_id_;
                continue;
            }
            props_data_.push_back(std::move(props_map));
            prop_names_.push_back(prop_itr.name());
        }
    }

    using view_data_t = std::map<double, std::vector<const Json::Value*>>;
    view_data_t view_data_;
    std::vector<std::map<double, Json::Value>> props_data_;
    std::vector<std::string> prop_names_;
    const Json::Value empty_val_;
    std::string layer_id_;
    double minzoom_;
    double maxzoom_;
    view_data_t::iterator iter_;
    bool valid_{false};
};

static bool ParseStyle(mapnik::Map &map, const std::string &style_str, const std::string& base_path);

class PropertyParser;

struct PatternLineParser;
struct LineParser;
struct BuildingParser;
struct PatternPolygonParser;
struct PolygonParser;
struct SymbolParser;
static std::string NormalizeFontName(const std::string& font_name);
static mapnik::formatting::node_ptr ParseTextField(const std::string& text_field);
static mapnik::expression_ptr ParseFilter(const Json::Value& mapbox_filter);


namespace detail {

template <typename T>
bool validate_property(const Json::Value& jprop) {
    LOG(ERROR) << "Unsupported requested property type " << typeid(T).name();
    return false;
}

template <>
inline bool validate_property<std::string>(const Json::Value& jprop) {
    return jprop.isString();
}

template <>
inline bool validate_property<bool>(const Json::Value& jprop) {
    return jprop.isBool();
}

template <>
inline bool validate_property<double>(const Json::Value& jprop) {
    return jprop.isDouble();
}

template <>
inline bool validate_property<mapnik::value_integer>(const Json::Value& jprop) {
    return jprop.isIntegral();
}

template <typename T>
inline T convert_property(const Json::Value& jprop) {
    return T();
}

template <>
inline std::string convert_property<std::string>(const Json::Value& jprop) {
    return jprop.asString();
}

template <>
inline bool convert_property<bool>(const Json::Value& jprop) {
    return jprop.asBool();
}

template <>
inline double convert_property<double>(const Json::Value& jprop) {
    return jprop.asDouble();
}

template <>
inline mapnik::value_integer convert_property<mapnik::value_integer>(const Json::Value& jprop) {
    return jprop.asInt64();
}

template <typename T>
boost::optional<T> get_property_impl(const std::string& prop_name, const Json::Value& jprop) {
    if (!validate_property<T>(jprop)) {
        LOG(ERROR) << "Property " << prop_name <<  ": " << jprop << " is not convertible to " << typeid(T).name();
        return boost::none;
    }
    return convert_property<T>(jprop);
}

template <>
boost::optional<mapnik::color> get_property_impl<mapnik::color>(const std::string& prop_name,
                                                                const Json::Value& jprop) {
    if (!jprop.isString()) {
        LOG(ERROR) << "Property " << prop_name <<  ": " << jprop << " must be a string!";
        return boost::none;
    }
    try {
        return mapnik::color(jprop.asString());
    } catch (const std::exception& e) {
        LOG(ERROR) << "Unable to parse color: " << jprop << ": " << e.what();
        return boost::none;
    }
}

template <>
boost::optional<mapnik::dash_array> get_property_impl<mapnik::dash_array>(const std::string& prop_name,
                                                                          const Json::Value& jprop) {
    if (jprop.isArray()) {
        Json::ArrayIndex len = jprop.size();
        if (len > 0 && len % 2 == 0) {
            mapnik::dash_array darray;
            Json::ArrayIndex i = 0;
            for (; i < len; i+=2) {
                const Json::Value& first = jprop[i];
                const Json::Value& second = jprop[i+1];
                if (!(first.isDouble() && second.isDouble())) {
                    break;
                }
                darray.emplace_back(first.asDouble(), second.asDouble());
            }
            // If all values are double
            if (i == len) {
                return darray;
            }
        }
    }
    LOG(ERROR) << "Property " << prop_name <<  ": " << jprop
               << " must be an array of double with size divisible by two!";
    return boost::none;
}

template <typename T>
boost::optional<T> get_property(const PropertiesView& prop_view, const std::string& prop_name) {
    const Json::Value& jprop = prop_view.property(prop_name);
    if (jprop.isNull()) {
        return boost::none;
    }
    return get_property_impl<T>(prop_name, jprop);
}

template <typename Variant>
Variant expr_node_to_variant(const mapnik::expr_node& expr) {
    return expr;
}

template <>
mapnik::symbolizer_base::value_type expr_node_to_variant(const mapnik::expr_node& expr) {
    return std::make_shared<mapnik::expr_node>(expr);
}

template <typename Variant, typename T>
boost::optional<Variant> get_variant(const PropertiesView& prop_view, const std::string& prop_name) {
    const Json::Value& jprop = prop_view.property(prop_name);
    if (jprop.isNull()) {
        return boost::none;
    }
    if (jprop.isObject()) {
        const Json::Value& jattr = jprop["property"];
        if (!jattr.isString()) {
            LOG(ERROR) << "Invalid property" << prop_name << ": " << jprop;
            return boost::none;
        }
        const Json::Value& jtype = jprop["type"];
        if (!jtype.isNull() && jtype != "identity") {
            LOG(ERROR) << "Unsupported function type " << jtype << " of property " << prop_name;
            return boost::none;
        }
        const std::string attr_name = jattr.asString();
        if (attr_name.empty()) {
            LOG(ERROR) << "Empty attribute name of property" << prop_name;
            return boost::none;
        }
        return expr_node_to_variant<Variant>(mapnik::attribute(attr_name));
    }

    auto val = get_property_impl<T>(prop_name, jprop);
    if (val) {
        return Variant(*val);
    }
    return boost::none;
}

} // ns detail


class PropertyParser {
public:

    using value_type = mapnik::symbolizer_base::value_type;

    PropertyParser(const PropertiesView& prop_view) : prop_view_(prop_view) {}

    template <typename T>
    boost::optional<T> get_property(const std::string& prop_name) const {
        return detail::get_property<T>(prop_view_, prop_name);
    }

    template <typename T>
    boost::optional<value_type> get_value(const std::string& prop_name) const {
        return detail::get_variant<value_type, T>(prop_view_, prop_name);
    }

    template <typename T>
    boost::optional<mapnik::expr_node> get_expr_node(const std::string& prop_name) const {
        return detail::get_variant<mapnik::expr_node, T>(prop_view_, prop_name);
    }

    template <typename T>
    void put_to_symbolizer(mapnik::symbolizer_base& sym, mapnik::keys key, const std::string& prop_name) const {
        const auto val = get_value<T>(prop_name);
        if (val) {
            mapnik::put(sym, key, *val);
        }
    }

    template <typename T>
    void assign_to_value(const std::string& prop_name, value_type& output_val) const {
        const auto val = get_value<T>(prop_name);
        if (val) {
            output_val = *val;
        }
    }

    boost::optional<std::string> get_filepath(const std::string& prop_name, const std::string& base_path) const {
        const auto file_path = get_property<std::string>(prop_name);
        if (file_path) {
            const std::string pattern_full_path = base_path + *file_path;
            if (boost::filesystem::exists(pattern_full_path)) {
                return pattern_full_path;
            } else {
                LOG(ERROR) << "File " << pattern_full_path << " not found!";
                return boost::none;
            }
        }
        return boost::none;
    }

    bool has_property(const std::string& prop_name) const {
        return !prop_view_.property(prop_name).isNull();
    }

    template <typename T>
    boost::optional<std::vector<T>> get_vector(const std::string& prop_name) const {
        const Json::Value& jprop = prop_view_.property(prop_name);
        if (jprop.isNull()) {
            return boost::none;
        }
        if (jprop.isArray()) {
            std::vector<T> result;
            Json::ArrayIndex len = jprop.size();
            result.reserve(len);
            bool success = true;
            for (Json::ArrayIndex i = 0; i < len; ++i) {
                const Json::Value& jitem = jprop[i];
                if (!detail::validate_property<T>(jitem)) {
                    success = false;
                    break;
                }
                result.push_back(detail::convert_property<T>(jitem));
            }
            if (success) {
                return result;
            }
        }
        LOG(ERROR) << "Property " << prop_name << ": " << jprop
                   << " must be an array of values convertible to " << typeid(T).name();
        return boost::none;
    }

private:
    const PropertiesView& prop_view_;
};

template <typename LayerParser>
void parse_layer(const Json::Value& mapbox_layer, const std::string& id,
                 std::vector<mapnik::rule> &rules, const std::string &res_path) {
    PropertiesView prop_view(mapbox_layer);
    while (prop_view.next()) {
        mapnik::rule r;
        r.set_name(id);
        r.set_min_scale(mapbox_zoom2scale(prop_view.maxzoom()));
        r.set_max_scale(mapbox_zoom2scale(prop_view.minzoom()));
        LayerParser parser;
        parser(PropertyParser(prop_view), res_path, r);
        rules.emplace_back(std::move(r));
    }
}

void load_mapbox_map(mapnik::Map & map, std::string const& filename, bool strict, const std::string& base_path) {
    std::ifstream file_stream(filename, std::ifstream::binary);
    if (!file_stream.is_open()) {
        throw std::runtime_error("Error openning file: " + filename);
    }
    std::string style_str(std::istreambuf_iterator<char>{file_stream}, {});
    const std::string& dir_path = base_path.empty() ? mapnik::util::dirname(filename) : base_path;
    load_mapbox_map_string(map, style_str, strict, dir_path);
}

void load_mapbox_map_string(mapnik::Map & map, std::string const& str, bool strict, const std::string& base_path) {
    if (!ParseStyle(map, str, base_path)) {
        throw std::runtime_error("Error while loading map!");
    }
}

static inline void WarnNotSupported(const PropertyParser& parser, const std::string& property_name) {
    if (parser.has_property(property_name)) {
        LOG(WARNING) << property_name << " not supported yet!";
    }
}


static bool ParseStyle(mapnik::Map &map, const std::string& style_str, const std::string& base_path) {
    Json::Value root;
    Json::CharReaderBuilder cr_builder;
    std::unique_ptr<Json::CharReader> char_reader(cr_builder.newCharReader());
    std::string json_erros;
    bool parsed = char_reader->parse(style_str.data(), style_str.data() + style_str.size(),
                                                &root, &json_erros);
    if (!parsed) {
        LOG(ERROR) << "Errors while pasing style json:\n" << json_erros;
        return false;
    }
    char_reader.reset();

    const Json::Value& root_layers = root["layers"];
    if (!root_layers.isArray()) {
        return false;
    }

    const std::string fonts_full_path = mapnik::util::make_absolute("fonts/", base_path);
    const std::string res_full_path = mapnik::util::make_absolute("res/", base_path);

    std::vector<mapnik::layer> mapnik_layers;
    // Map in which key is the name of the source layer and
    // value is the style vector associated with this source layer
    std::map<std::string, std::vector<mapnik::feature_type_style>> styles_map;

    for (unsigned int layer_num = 0; layer_num < root_layers.size(); ++layer_num) {
        const Json::Value& mapbox_layer = root_layers[layer_num];
        if (!mapbox_layer.isObject()) {
            LOG(ERROR) << "Invalid layer!\n" << mapbox_layer;
            continue;
        }
        const Json::Value& mapbox_paint = mapbox_layer["paint"];
        std::string id = mapbox_layer.get("id", "").asString();
        std::string source_layer_name = mapbox_layer.get("source-layer", "").asString();
        std::string layer_type = mapbox_layer.get("type", "").asString();

        std::vector<mapnik::rule> rules;

        if (layer_type == "background") {
            std::string bgcolor_string = mapbox_paint.get("background-color", "#000000").asString();
            mapnik::color bgcolor(bgcolor_string);
            map.set_background(bgcolor);
            continue;
        } else if (layer_type == "fill") {
            if (mapbox_paint.isMember("fill-pattern")) {
                parse_layer<PatternPolygonParser>(mapbox_layer, id, rules, res_full_path);
            } else {
                parse_layer<PolygonParser>(mapbox_layer, id, rules, res_full_path);
            }
        } else if (layer_type == "fill-extrusion") {
            parse_layer<BuildingParser>(mapbox_layer, id, rules, res_full_path);
        } else if (layer_type == "line") {
            if (mapbox_paint.isMember("line-pattern")) {
                parse_layer<PatternLineParser>(mapbox_layer, id, rules, res_full_path);
            } else {
                parse_layer<LineParser>(mapbox_layer, id, rules, res_full_path);
            }
        } else if (layer_type == "symbol") {
            parse_layer<SymbolParser>(mapbox_layer, id, rules, res_full_path);
        } else {
            LOG(ERROR) << "Invalid layer type: " << layer_type;
        }

        if (rules.empty()) {
            continue;
        }

        double mapbox_minzoom = mapbox_layer.get("minzoom", 0).asDouble();
        double mapbox_maxzoom = mapbox_layer.get("maxzoom", 20).asDouble();
        double current_min_sd = mapbox_zoom2scale(mapbox_maxzoom);
        double current_max_sd = mapbox_zoom2scale(mapbox_minzoom);

        auto styles_itr = styles_map.find(source_layer_name);
        if (styles_itr == styles_map.end()) {
            // No styles associated with this source layer,
            // so we should create new vector of styles containing one style
            auto emplace_res = styles_map.emplace(source_layer_name, std::vector<mapnik::feature_type_style>());
            styles_itr = emplace_res.first;
        }

        auto& styles_vector = styles_itr->second;

        if (!mapnik_layers.empty() && mapnik_layers.back().name() == source_layer_name) {
            // Current and previous mapbox layers has the same source layer, so we should union them into one layer
            mapnik::layer& current_layer = mapnik_layers.back();
            current_layer.set_maximum_scale_denominator(std::max(current_layer.maximum_scale_denominator(),
                                                                 current_max_sd));
            current_layer.set_minimum_scale_denominator(std::min(current_layer.minimum_scale_denominator(),
                                                                 current_min_sd));
        } else {
            // Here we should create new layer and style
            mapnik_layers.emplace_back(source_layer_name, mapnik::MAPNIK_GMERC_PROJ);
            mapnik::layer& new_layer = mapnik_layers.back();
            new_layer.set_minimum_scale_denominator(current_min_sd);
            new_layer.set_maximum_scale_denominator(current_max_sd);
        }

        // Append new style to styles vector
        styles_vector.emplace_back();
        uint style_num = styles_vector.size() - 1;
        mapnik_layers.back().add_style(source_layer_name + std::to_string(style_num));

        // Current style should be the last in current styles vector
        mapnik::feature_type_style& current_style = styles_vector.back();
        const Json::Value& mapbox_filter = mapbox_layer["filter"];
        for (auto& rule : rules) {
            if (!mapbox_filter.isNull()) {
                mapnik::expression_ptr f = ParseFilter(mapbox_filter);
                rule.set_filter(f);
            }
            current_style.add_rule(std::move(rule));
        }
    }

    map.set_font_directory(fonts_full_path);
    map.register_fonts(fonts_full_path);

    //Insert all layers into map
    for (auto& lyr : mapnik_layers) {
        map.add_layer(std::move(lyr));
    }

    // Insert all styles into map
    for (auto styles_map_itr = styles_map.begin(); styles_map_itr != styles_map.end(); ++styles_map_itr)
    {
        const std::string& source_layer_name = styles_map_itr->first;
        auto& styles_vector = styles_map_itr->second;
        auto num_styles = styles_vector.size();
        for (uint i = 0; i < num_styles; ++i) {
            map.insert_style(source_layer_name + std::to_string(i), std::move(styles_vector[i]));
        }
    }

    return true;
}

using keys = mapnik::keys;
using vertical_alignment_enum = mapnik::vertical_alignment_enum;
using horizontal_alignment_enum = mapnik::horizontal_alignment_enum;
using justify_alignment_enum = mapnik::justify_alignment_enum;

struct PatternLineParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {
        const auto pattern_full_path = prop_parser.get_filepath("line-pattern", res_path);
        if (!pattern_full_path) {
            LOG(ERROR) << "No pattern found!";
            return;
        }
        mapnik::line_pattern_symbolizer lps;
        mapnik::put(lps, mapnik::keys::file, *pattern_full_path);
        output_rule.append(std::move(lps));
    }
};

struct LineParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {
        mapnik::line_symbolizer ls;

        prop_parser.put_to_symbolizer<mapnik::color>(ls, keys::stroke, "line-color");
        prop_parser.put_to_symbolizer<double>(ls, keys::stroke_opacity, "line-opacity");

        double line_width = 1.0;
        const auto lwidth = prop_parser.get_property<double>("line-width");
        if (lwidth) {
            line_width = *lwidth;
            mapnik::put(ls, keys::stroke_width, *lwidth);
        }

        const auto ljoin = prop_parser.get_property<std::string>("line-join");
        if (ljoin) {
            const std::string& join_type = *ljoin;
            std::uint8_t mapnik_join_type = mapnik::MITER_JOIN;
            if (join_type == "round") {
                mapnik_join_type = mapnik::ROUND_JOIN;
            } else if(join_type == "miter_revert") {
                mapnik_join_type= mapnik::MITER_REVERT_JOIN;
            } else if(join_type == "bevel") {
                mapnik_join_type = mapnik::BEVEL_JOIN;
            } else {
                LOG(ERROR) << "Invalid join type: '" << join_type << "'!";
            }
            mapnik::put(ls, keys::stroke_linejoin, mapnik::line_join_enum(mapnik_join_type));
        }

        const auto lcap = prop_parser.get_property<std::string>("line-cap");
        if (lcap) {
            const std::string& join_type = *lcap;
            std::uint8_t mapnik_cap_type = mapnik::BUTT_CAP;
            if (join_type == "round") {
                mapnik_cap_type = mapnik::ROUND_CAP;
            } else if(join_type == "square") {
                mapnik_cap_type= mapnik::SQUARE_CAP;
            } else if(join_type == "butt") {
                mapnik_cap_type= mapnik::BUTT_CAP;
            } else {
                LOG(ERROR) << "Invalid join type: '" << join_type << "'!";
            }
            mapnik::put(ls, keys::stroke_linecap, mapnik::line_cap_enum(mapnik_cap_type));
        }

        auto ldash = prop_parser.get_property<mapnik::dash_array>("line-dasharray");
        if (ldash) {
            for (auto& p : *ldash) {
                p.first *= line_width;
                p.second *= line_width;
            }
            mapnik::put(ls, keys::stroke_dasharray, *ldash);
        }

        // Should be the last property checked
        const auto lgwidth = prop_parser.get_property<double>("line-gap-width");
        if (lgwidth) {
            mapnik::line_symbolizer other_side(ls);
            double half_lgwidth = *lgwidth / 2;
            mapnik::put(ls, keys::offset, half_lgwidth);
            mapnik::put(other_side, mapnik::keys::offset, -half_lgwidth);
            output_rule.append(std::move(other_side));
        }

        output_rule.append(std::move(ls));
    }
};

struct BuildingParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {
        mapnik::building_symbolizer bs;
        prop_parser.put_to_symbolizer<double>(bs, mapnik::keys::height, "fill-extrusion-heigh");
        prop_parser.put_to_symbolizer<mapnik::color>(bs, mapnik::keys::fill, "fill-extrusion-color");
        prop_parser.put_to_symbolizer<double>(bs, mapnik::keys::fill_opacity, "fill-extrusion-opacity");
        output_rule.append(std::move(bs));
    }
};

struct PatternPolygonParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {
        const auto pattern_full_path = prop_parser.get_filepath("fill-pattern", res_path);
        if (!pattern_full_path) {
            LOG(ERROR) << "No pattern found!";
            return;
        }
        mapnik::polygon_pattern_symbolizer pps;
        mapnik::put(pps, mapnik::keys::file, *pattern_full_path);
        prop_parser.put_to_symbolizer<double>(pps, mapnik::keys::fill_opacity, "fill-opacity");
        output_rule.append(std::move(pps));
    }
};

struct PolygonParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {
        mapnik::polygon_symbolizer ps;
        prop_parser.put_to_symbolizer<mapnik::color>(ps, mapnik::keys::fill, "fill-color");
        prop_parser.put_to_symbolizer<double>(ps, mapnik::keys::fill_opacity, "fill-opacity");
        output_rule.append(std::move(ps));
    }
};


static void ParseTextProperties(const PropertyParser &prop_parser, mapnik::symbolizer_base& symb) {
    double text_size = kEmsInPixel;
    const auto tsize = prop_parser.get_property<double>("text-size");
    if (tsize) {
        text_size = *tsize;
    }

    mapnik::text_placements_ptr p;
    bool use_simple_placement = false;

    const auto placements = prop_parser.get_property<std::string>("me:placements");
#ifdef MAPBOX2MAPNIK_ENABLE_SIMPLE_PLACEMENTS
    if (placements) {
        use_simple_placement = true;
        // TODO: validate string
        p = std::make_shared<mapnik::text_placements_simple>(*placements);

        const auto dx = prop_parser.get_property<double>("me:dx");
        if (dx) {
            p->defaults.layout_defaults.dx = ems2pixels(*dx, text_size);
        }
        const auto dy = prop_parser.get_property<double>("me:dy");
        if (dy) {
            p->defaults.layout_defaults.dy = ems2pixels(*dy, text_size);
        }
    } else {
        p = std::make_shared<mapnik::text_placements_dummy>();
    }
#else
    if (placements) {
        LOG(ERROR) << "Simple placements not supported in this build!";
    }
    p = std::make_shared<mapnik::text_placements_dummy>();
#endif

    p->defaults.format_defaults.text_size = text_size;

    p->defaults.expressions.avoid_edges = true;
    p->defaults.expressions.allow_overlap = false;
    //should be done because not everywhere font is determined
    //p->defaults.format_defaults.face_name = "DejaVu Sans Bold";

    // the text to be displayed
    auto tfield = prop_parser.get_property<std::string>("text-field");
    if (tfield) {
        p->defaults.set_format_tree(ParseTextField(*tfield));
    }

    const auto tfonts = prop_parser.get_vector<std::string>("text-font");
    if (tfonts && (*tfonts).size() > 0) {
        p->defaults.format_defaults.face_name = NormalizeFontName((*tfonts)[0]);
    }

    prop_parser.assign_to_value<double>("text-letter-spacing", p->defaults.format_defaults.character_spacing);
    prop_parser.assign_to_value<mapnik::color>("text-color", p->defaults.format_defaults.fill);
    prop_parser.assign_to_value<double>("text-opacity", p->defaults.format_defaults.text_opacity);
    prop_parser.assign_to_value<mapnik::color>("text-halo-color", p->defaults.format_defaults.halo_fill);
    prop_parser.assign_to_value<double>("text-halo-width", p->defaults.format_defaults.halo_radius);
    prop_parser.assign_to_value<double>("symbol-spacing", p->defaults.expressions.repeat_distance);
    prop_parser.assign_to_value<double>("text-max-angle", p->defaults.expressions.max_char_angle_delta);

    const auto placement = prop_parser.get_property<std::string>("symbol-placement");
    if (placement) {
        const std::string& placement_str = *placement;
        std::uint8_t marker_placement = mapnik::POINT_PLACEMENT;
        if (placement_str == "point") {
            marker_placement = mapnik::POINT_PLACEMENT;
        } else if (placement_str == "line") {
            marker_placement = mapnik::LINE_PLACEMENT;
        } else if(placement_str == "vertex") {
            marker_placement = mapnik::VERTEX_PLACEMENT;
        } else if(placement_str == "interior") {
            marker_placement = mapnik::INTERIOR_PLACEMENT;
        } else {
            LOG(ERROR) << "Invalid placement type: '" << placement_str << "'!";
        }
        p->defaults.expressions.label_placement = mapnik::enumeration_wrapper(marker_placement);
    }

    const auto tpading = prop_parser.get_property<double>("text-padding");
    const auto spading = prop_parser.get_property<double>("symbol-padding");
    if(tpading || spading) {
        if(tpading) {
            if (spading) {
                p->defaults.expressions.margin = std::max(*tpading, *spading);
            } else {
                p->defaults.expressions.margin = *tpading;
            }
        } else {
            p->defaults.expressions.margin = *spading;
        }
    };

    const auto tmwidth = prop_parser.get_property<double>("text-max-width");
    if (tmwidth) {
        p->defaults.layout_defaults.wrap_before = true;
        p->defaults.layout_defaults.wrap_width = ems2pixels(*tmwidth, text_size);
    }

    const auto tlheight = prop_parser.get_property<double>("text-line-height");
    if (tlheight) {
        p->defaults.format_defaults.line_spacing = ems2pixels(*tlheight, text_size);
    }

    const auto tjustify = prop_parser.get_property<std::string>("text-justify");
    if (tjustify) {
        const std::string& justify_str = *tjustify;
        std::uint8_t jalign = mapnik::J_AUTO;
        if (justify_str == "left") {
            jalign = mapnik::J_LEFT;
        } else if (justify_str == "right") {
            jalign = mapnik::J_RIGHT;
        } else if (justify_str == "center") {
            jalign = mapnik::J_MIDDLE;
        } else {
            LOG(ERROR) << "Invalid justify type: '" << justify_str << "'!";
        }
        p->defaults.layout_defaults.jalign = mapnik::enumeration_wrapper(jalign);
    }

    if (!use_simple_placement) {
        const auto toffset = prop_parser.get_vector<double>("text-offset");
        if (toffset) {
            const auto offset_vec = *toffset;
            std::size_t len = offset_vec.size();
            if (len > 0) {
                p->defaults.layout_defaults.dx = ems2pixels(offset_vec[0], text_size);
            }
            if (len > 1) {
                p->defaults.layout_defaults.dy = ems2pixels(offset_vec[1], text_size);
            }
        }

        std::uint8_t v_align = mapnik::V_MIDDLE;
        std::uint8_t h_align = mapnik::V_MIDDLE;
        const auto tanchor = prop_parser.get_property<std::string>("text-anchor");
        if (tanchor) {
            const std::string& anchor_str = *tanchor;
            if (anchor_str.find("center") == anchor_str.npos) {
                if (anchor_str.find("top") != anchor_str.npos) {
                    v_align = mapnik::V_BOTTOM;
                } else if (anchor_str.find("bottom") != anchor_str.npos) {
                    v_align = mapnik::V_TOP;
                }
                if (anchor_str.find("left") != anchor_str.npos) {
                    h_align = mapnik::H_RIGHT;
                } else if (anchor_str.find("right") != anchor_str.npos) {
                    h_align = mapnik::H_LEFT;
                }
            }

        }
        p->defaults.layout_defaults.valign = mapnik::enumeration_wrapper(v_align);
        p->defaults.layout_defaults.halign = mapnik::enumeration_wrapper(h_align);
    }

    WarnNotSupported(prop_parser, "text-transform");
    WarnNotSupported(prop_parser, "text-ignore-placement");
    WarnNotSupported(prop_parser, "text-halo-blur");

    mapnik::put(symb, keys::text_placements_, p);
}

static void ParseIconProperties(const PropertyParser &prop_parser, const std::string& res_path,
                         mapnik::symbolizer_base& symb) {
    const auto full_icon_path = prop_parser.get_filepath("icon-image", res_path);
    if (full_icon_path) {
        mapnik::put(symb, mapnik::keys::file, *full_icon_path);
    }

    mapnik::transform_list_ptr t_list_ptr = std::make_shared<mapnik::transform_list>();
    const auto isize = prop_parser.get_expr_node<double>("icon-size");
    if (isize) {
        t_list_ptr->push_back(mapnik::scale_node(*isize, boost::none));
    }

    const auto irotate = prop_parser.get_expr_node<double>("icon-rotate");
    if (irotate) {
        t_list_ptr->push_back(mapnik::rotate_node(*irotate));
    }

    if (!t_list_ptr->empty()) {
        mapnik::put(symb, mapnik::keys::image_transform, t_list_ptr);
    }

    prop_parser.put_to_symbolizer<double>(symb, mapnik::keys::opacity, "icon-opacity");

    WarnNotSupported(prop_parser, "icon-text-fit");
}


static void ParseMarkersSymb(const PropertyParser &prop_parser, const std::string& res_path,
                             mapnik::symbolizer_base& symb) {
    ParseIconProperties(prop_parser, res_path, symb);
    prop_parser.put_to_symbolizer<double>(symb, mapnik::keys::spacing, "symbol-spacing");
    mapnik::put(symb, mapnik::keys::avoid_edges, true);

    const auto placement = prop_parser.get_property<std::string>("symbol-placement");
    if (placement) {
        const std::string& placement_str = *placement;
        std::uint8_t marker_placement = mapnik::MARKER_POINT_PLACEMENT;
        if (placement_str == "point") {
            marker_placement = mapnik::MARKER_POINT_PLACEMENT;
        } else if (placement_str == "line") {
            marker_placement = mapnik::MARKER_LINE_PLACEMENT;
        } else if(placement_str == "vertex") {
                marker_placement = mapnik::MARKER_VERTEX_FIRST_PLACEMENT;
        } else if(placement_str == "interior") {
                marker_placement = mapnik::MARKER_INTERIOR_PLACEMENT;
        } else {
            LOG(ERROR) << "Invalid placement type: '" << placement_str << "'!";
        }
        put(symb, mapnik::keys::markers_placement_type, mapnik::marker_placement_enum(marker_placement));
    }
}

static void ParseShieldSymb(const PropertyParser &prop_parser, const std::string& res_path,
                             mapnik::symbolizer_base& symb) {
    ParseTextProperties(prop_parser, symb);
    ParseIconProperties(prop_parser, res_path, symb);
    mapnik::put(symb, mapnik::keys::unlock_image, true);
}

static void ParseTextSymb(const PropertyParser &prop_parser, mapnik::symbolizer_base& symb) {
    ParseTextProperties(prop_parser, symb);
}

struct SymbolParser {
    void operator() (const PropertyParser &prop_parser, const std::string &res_path, mapnik::rule& output_rule) {

        bool has_text = prop_parser.has_property("text-field");
        bool has_image = prop_parser.has_property("icon-image");

        mapnik::symbolizer_base symb;

        if (has_image) {
            if (has_text) {
                ParseShieldSymb(prop_parser, res_path, symb);
                output_rule.append(static_cast<mapnik::shield_symbolizer&&>(std::move(symb)));
            } else {
                ParseMarkersSymb(prop_parser, res_path, symb);
                output_rule.append(static_cast<mapnik::markers_symbolizer&&>(std::move(symb)));
            }
        } else {
            if (has_text) {
                ParseTextSymb(prop_parser, symb);
                output_rule.append(static_cast<mapnik::text_symbolizer&&>(std::move(symb)));
            } else {
                LOG(ERROR) << "Symbol does not have either a text or an image!";
                return;
            }
        }
    }
};


using or_node = mapnik::binary_node<mapnik::tags::logical_or>;
using and_node = mapnik::binary_node<mapnik::tags::logical_and>;
using equal_to_node = mapnik::binary_node<mapnik::tags::equal_to>;
using not_equal_to_node = mapnik::binary_node<mapnik::tags::not_equal_to>;

static mapnik::expr_node make_filter_node(const Json::Value& mapbox_filter);

template <typename LogicalNode>
static inline mapnik::expr_node combine_expr_nodes(const Json::ValueConstIterator& begin,
                                                   const Json::ValueConstIterator& end) {
    auto iter = begin;
    mapnik::expr_node result = make_filter_node(*iter);
    for (++iter; iter != end; ++iter) {
        result = LogicalNode(result, make_filter_node(*iter));
    }
    return result;
}

static mapnik::expr_node make_filter_node(const Json::Value& mapbox_filter) {
    auto val_type = mapbox_filter.type();
    if (val_type == Json::booleanValue) {
        return mapbox_filter.asBool();
    }
    if (val_type == Json::intValue) {
        return static_cast<mapnik::value_integer>(mapbox_filter.asInt());
    }
    if (val_type == Json::realValue) {
        return mapbox_filter.asDouble();
    }
    if (val_type == Json::stringValue) {
        const std::string val = mapbox_filter.asString();
        return str2expr_node(val);
    }
    if (mapbox_filter.isUInt64()) {
        return static_cast<mapnik::value_integer>(mapbox_filter.asUInt64());
    }

    uint filter_size = mapbox_filter.size();
    if (!mapbox_filter.isArray() && filter_size < 2) {
        LOG(ERROR) << "Invalid filter!" << std::endl << mapbox_filter;
        return mapnik::value_null();
    }

    auto first_val_itr = mapbox_filter.begin();
    std::string filter_type = first_val_itr->asString();
    ++first_val_itr;
    if (filter_type == "all") {
        return combine_expr_nodes<and_node>(first_val_itr, mapbox_filter.end());
    }
    if (filter_type == "any") {
        return combine_expr_nodes<or_node>(first_val_itr, mapbox_filter.end());
    }
    if (filter_type == "none") {
        return  mapnik::unary_node<mapnik::tags::logical_not>(
                    combine_expr_nodes<and_node>(first_val_itr, mapbox_filter.end()));
    }

    // Attribute dependent filters
    const Json::Value& jattr = mapbox_filter[1];
    if (!jattr.isString()) {
        LOG(ERROR) << "Invalid attribute " << jattr << " in filter " << mapbox_filter;
        return mapnik::value_null();
    }
    mapnik::attribute attr(jattr.asString());
    if (filter_type == "has") {
        return not_equal_to_node(attr, mapnik::value_null());
    }
    if (filter_type == "!has") {
        return equal_to_node(attr, mapnik::value_null());
    }
    if (filter_size > 2) {
        if (filter_type == "in") {
            mapnik::expr_node result = equal_to_node(attr, make_filter_node(mapbox_filter[2]));
            for (uint i = 3; i < mapbox_filter.size(); ++i) {
                result = or_node(result, equal_to_node(attr, make_filter_node(mapbox_filter[i])));
            }
            return result;
        }
        if (filter_type == "!in") {
            mapnik::expr_node result = not_equal_to_node(attr, make_filter_node(mapbox_filter[2]));
            for (uint i = 3; i < mapbox_filter.size(); ++i) {
                result = and_node(result, not_equal_to_node(attr, make_filter_node(mapbox_filter[i])));
            }
            return result;
        }
        if (filter_type == "==") {
            return equal_to_node(attr, make_filter_node(mapbox_filter[2]));
        }
        if (filter_type == "!=") {
            return mapnik::binary_node<mapnik::tags::not_equal_to>(attr, make_filter_node(mapbox_filter[2]));
        }
        if (filter_type == ">") {
            return mapnik::binary_node<mapnik::tags::greater>(attr, make_filter_node(mapbox_filter[2]));
        }
        if (filter_type == ">=") {
            return mapnik::binary_node<mapnik::tags::greater_equal>(attr, make_filter_node(mapbox_filter[2]));
        }
        if (filter_type == "<") {
            return mapnik::binary_node<mapnik::tags::less>(attr, make_filter_node(mapbox_filter[2]));
        }
        if (filter_type == "<=") {
            return mapnik::binary_node<mapnik::tags::less_equal>(attr, make_filter_node(mapbox_filter[2]));
        }
    }

    LOG(ERROR) << "Invalid filter!" << std::endl << mapbox_filter;
    return mapnik::value_null();
}

static mapnik::expression_ptr ParseFilter(const Json::Value& mapbox_filter) {
    if (!mapbox_filter.isArray() || mapbox_filter.size() < 2) {
        LOG(ERROR) << "Invalid filter!" << std::endl << mapbox_filter;
        return nullptr;
    }
    return std::make_shared<mapnik::expr_node>(make_filter_node(mapbox_filter));
}


static std::string NormalizeFontName(const std::string& font_name) {
    // TODO: make generic
    if(font_name == "DejaVu Sans Bold" || font_name == "DejaVuSans-Bold")
        return "DejaVu Sans Bold";
    else if(font_name == "DejaVuSans")
        return "DejaVu Sans Book";
    else if(font_name == "DejaVuSansCondensed")
        return "DejaVu Sans Condensed";
    else if (font_name == "DejaVuSans-Oblique")
        return "DejaVu Sans Oblique";
    else if (font_name == "DejaVuSansCondensed-Oblique")
        return "DejaVu Sans Condensed Oblique";
    return "DejaVu Sans Bold";
}

static mapnik::formatting::node_ptr ParseTextField(const std::string& text_field) {
    std::size_t tf_len = text_field.size();
    if (tf_len > 1) {
        std::size_t s_pos = text_field.find('{');
        if (s_pos != text_field.npos) {
            std::size_t expr_pos = s_pos + 1;
            std::size_t e_pos = text_field.find('}', expr_pos);
            if (e_pos != text_field.npos) {
                using plus_node = mapnik::binary_node<mapnik::tags::plus>;
                mapnik::expr_node expr;
                expr = mapnik::attribute(text_field.substr(expr_pos, e_pos - expr_pos));
                if (s_pos > 0) {
                    expr = plus_node(str2expr_node(text_field.substr(0, s_pos)), expr);
                }
                std::size_t search_start = e_pos + 1;
                while (search_start != tf_len) {
                    s_pos = text_field.find('{', search_start);
                    if (s_pos == text_field.npos) {
                        break;
                    }
                    if (s_pos > search_start) {
                        expr = plus_node(expr, str2expr_node(text_field.substr(search_start, s_pos - search_start)));
                    }
                    expr_pos = s_pos + 1;
                    e_pos = text_field.find('}', expr_pos);
                    if (e_pos == text_field.npos) {
                        break;
                    }
                    expr = plus_node(expr, mapnik::attribute(text_field.substr(expr_pos, e_pos - expr_pos)));
                    search_start = e_pos + 1;
                }
                if (search_start != tf_len) {
                    expr = plus_node(expr, str2expr_node(text_field.substr(search_start)));
                }
                return std::make_shared<mapnik::formatting::text_node>(
                            std::make_shared<mapnik::expr_node>(std::move(expr)));
            }
        }
    }
    return std::make_shared<mapnik::formatting::text_node>(
                std::make_shared<mapnik::expr_node>(str2expr_node(text_field)));
}

} // ns sputnik
