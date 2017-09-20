#pragma once

#include <boost/optional.hpp>

#include <json/value.h>

#include <glog/logging.h>


namespace sputnik {

namespace detail {

template <typename T>
bool validate(const Json::Value& jvalue);

template <>
inline bool validate<std::string>(const Json::Value& jvalue) {
    return jvalue.isString();
}

template <>
inline bool validate<int>(const Json::Value& jvalue) {
    return jvalue.isInt();
}

template <>
inline bool validate<std::int64_t>(const Json::Value& jvalue) {
    return jvalue.isInt64();
}

template <>
inline bool validate<uint>(const Json::Value& jvalue) {
    return jvalue.isUInt();
}

template <>
inline bool validate<std::uint64_t>(const Json::Value& jvalue) {
    return jvalue.isUInt64();
}


template <>
inline bool validate<double>(const Json::Value& jvalue) {
    return jvalue.isDouble();
}

template <>
inline bool validate<bool>(const Json::Value& jvalue) {
    return jvalue.isBool();
}

template <typename T>
T from_json(const Json::Value& jvalue);

template <>
inline std::string from_json<std::string>(const Json::Value& jvalue) {
    return jvalue.asString();
}

template <>
inline int from_json<int>(const Json::Value& jvalue) {
    return jvalue.asInt();
}

template <>
inline std::int64_t from_json<std::int64_t>(const Json::Value& jvalue) {
    return jvalue.asInt64();
}

template <>
inline uint from_json<uint>(const Json::Value& jvalue) {
    return jvalue.asUInt();
}

template <>
inline std::uint64_t from_json<std::uint64_t>(const Json::Value& jvalue) {
    return jvalue.asUInt64();
}

template <>
inline double from_json<double>(const Json::Value& jvalue) {
    return jvalue.asDouble();
}

template <>
inline bool from_json<bool>(const Json::Value& jvalue) {
    return jvalue.asBool();
}

} // ns detail

using boost::optional;

template <typename T>
inline T FromJson(const Json::Value& jvalue, const T& default_value) {
    if (!detail::validate<T>(jvalue)) {
        return default_value;
    }
    return detail::from_json<T>(jvalue);
}

template <typename T>
inline optional<T> FromJson(const Json::Value& jvalue) {
    if (!detail::validate<T>(jvalue)) {
        return boost::none;
    }
    return detail::from_json<T>(jvalue);
}


template <typename T>
inline optional<T> FromJsonOrErr(const Json::Value& jvalue, const std::string& err_string) {
    auto val = FromJson<T>(jvalue);
    if (!val) {
        LOG(ERROR) << err_string;
    }
    return val;
}

} // ns sputnik
