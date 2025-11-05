#pragma once

#define HSML_C_API_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
  #if defined(HSML_C_API_BUILD)
    #define HSML_API __declspec(dllexport)
  #else
    #define HSML_API __declspec(dllimport)
  #endif
#else
  #define HSML_API __attribute__((visibility("default")))
#endif

// Returns 1 on success, 0 on failure
HSML_API int hsml_init(void);

// Shutdown framework; safe to call multiple times
HSML_API void hsml_shutdown(void);

// Returns 1 if initialized, 0 otherwise
HSML_API int hsml_is_initialized(void);

// Writes null-terminated version string into buffer; returns bytes written (excluding null).
// If buffer_size == 0 or buffer == NULL, returns required size (excluding null).
HSML_API int hsml_get_version(char* buffer, int buffer_size);

// Writes null-terminated build info (date time) into buffer; same contract as hsml_get_version
HSML_API int hsml_get_build_info(char* buffer, int buffer_size);

// Validate HSML file; returns 1 if valid, 0 otherwise
HSML_API int hsml_validate_file(const char* path);

#ifdef __cplusplus
}
#endif


