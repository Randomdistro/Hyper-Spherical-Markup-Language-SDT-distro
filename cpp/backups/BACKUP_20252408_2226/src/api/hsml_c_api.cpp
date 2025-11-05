#include "hsml/c/hsml_c.h"
#include "hsml/hsml.h"
#include <cstring>
#include <cstdio>

extern "C" {

HSML_API int hsml_init(void) {
    hsml::Framework::InitOptions opts{};
    return hsml::Framework::initialize(opts) ? 1 : 0;
}

HSML_API void hsml_shutdown(void) {
    hsml::Framework::shutdown();
}

HSML_API int hsml_is_initialized(void) {
    return hsml::Framework::is_initialized() ? 1 : 0;
}

HSML_API int hsml_get_version(char* buffer, int buffer_size) {
    const char* v = hsml::FrameworkInfo::version_string;
    const int len = static_cast<int>(std::strlen(v));
    if (!buffer || buffer_size <= 0) return len;
    const int n = (len < buffer_size - 1) ? len : (buffer_size - 1);
    std::memcpy(buffer, v, n);
    buffer[n] = '\0';
    return n;
}

HSML_API int hsml_get_build_info(char* buffer, int buffer_size) {
    // date + space + time
    const char* date = hsml::FrameworkInfo::build_date;
    const char* time = hsml::FrameworkInfo::build_time;
    char tmp[128]{};
    std::snprintf(tmp, sizeof(tmp), "%s %s", date ? date : "", time ? time : "");
    const int len = static_cast<int>(std::strlen(tmp));
    if (!buffer || buffer_size <= 0) return len;
    const int n = (len < buffer_size - 1) ? len : (buffer_size - 1);
    std::memcpy(buffer, tmp, n);
    buffer[n] = '\0';
    return n;
}

HSML_API int hsml_validate_file(const char* path) {
    if (!path || !*path) return 0;
    hsml::parsers::HSMLParser parser;
    parser.set_strict_mode(true);
    parser.enable_semantic_validation(true);
    auto bubble = parser.parse_file(path);
    return (bubble && !parser.has_errors()) ? 1 : 0;
}

}


