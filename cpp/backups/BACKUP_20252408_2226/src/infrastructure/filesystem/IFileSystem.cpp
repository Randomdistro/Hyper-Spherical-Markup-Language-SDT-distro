/** @file IFileSystem.h
 * @brief Cross-platform file system abstraction
 *
 * Clean Architecture: Infrastructure Layer - File System Abstraction
 * Provides unified file system operations across platforms.
 */

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <chrono>
#include <functional>
#include <unordered_map>

namespace hsml {
namespace infrastructure {

// File attributes
struct FileAttributes {
    bool exists{false};
    bool is_directory{false};
    bool is_file{false};
    bool is_readable{false};
    bool is_writable{false};
    bool is_executable{false};
    bool is_hidden{false};
    uint64_t size{0};
    std::chrono::system_clock::time_point creation_time;
    std::chrono::system_clock::time_point last_access_time;
    std::chrono::system_clock::time_point last_write_time;
};

// File open modes
enum class FileOpenMode {
    READ = 1 << 0,
    WRITE = 1 << 1,
    APPEND = 1 << 2,
    CREATE = 1 << 3,
    TRUNCATE = 1 << 4,
    BINARY = 1 << 5,
    TEXT = 1 << 6
};

inline FileOpenMode operator|(FileOpenMode a, FileOpenMode b) {
    return static_cast<FileOpenMode>(static_cast<int>(a) | static_cast<int>(b));
}

inline FileOpenMode operator&(FileOpenMode a, FileOpenMode b) {
    return static_cast<FileOpenMode>(static_cast<int>(a) & static_cast<int>(b));
}

// Directory entry
struct DirectoryEntry {
    std::string name;
    std::string full_path;
    FileAttributes attributes;
};

// File watch event types
enum class FileWatchEventType {
    CREATED,
    DELETED,
    MODIFIED,
    RENAMED
};

// File watch event
struct FileWatchEvent {
    FileWatchEventType type;
    std::string path;
    std::string old_path; // For renames
    std::chrono::system_clock::time_point timestamp;
};

// Path types
enum class PathType {
    ABSOLUTE,
    RELATIVE,
    HOME,
    APP_DATA,
    DOCUMENTS,
    DESKTOP,
    DOWNLOADS,
    PICTURES,
    MUSIC,
    VIDEOS,
    CACHE,
    TEMP,
    CONFIG,
    LOGS,
    SAVE_DATA
};

/**
 * @brief File handle for platform-specific file operations
 */
class IFileHandle {
public:
    virtual ~IFileHandle() = default;

    virtual bool is_open() const = 0;
    virtual void close() = 0;
    virtual uint64_t read(void* buffer, uint64_t size) = 0;
    virtual uint64_t write(const void* buffer, uint64_t size) = 0;
    virtual bool seek(int64_t offset, int origin) = 0;
    virtual uint64_t tell() const = 0;
    virtual uint64_t size() const = 0;
    virtual bool flush() = 0;
    virtual bool eof() const = 0;
};

/**
 * @brief Cross-platform file system interface
 *
 * This interface provides unified file system operations across platforms,
 * supporting both synchronous and asynchronous operations with proper
 * error handling and spherical coordinate data persistence.
 */
class IFileSystem {
public:
    virtual ~IFileSystem() = default;

    // Initialization
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Path operations
    virtual std::string get_current_directory() const = 0;
    virtual bool set_current_directory(const std::string& path) = 0;
    virtual std::string get_executable_path() const = 0;
    virtual std::string get_user_home_directory() const = 0;
    virtual std::string get_special_path(PathType type) const = 0;

    // Path manipulation
    virtual std::string normalize_path(const std::string& path) const = 0;
    virtual std::string absolute_path(const std::string& path) const = 0;
    virtual std::string relative_path(const std::string& from, const std::string& to) const = 0;
    virtual std::string combine_path(const std::string& base, const std::string& relative) const = 0;
    virtual std::string get_filename(const std::string& path) const = 0;
    virtual std::string get_directory(const std::string& path) const = 0;
    virtual std::string get_extension(const std::string& path) const = 0;
    virtual std::string change_extension(const std::string& path, const std::string& new_ext) const = 0;
    virtual bool has_extension(const std::string& path, const std::string& ext) const = 0;

    // File operations
    virtual bool file_exists(const std::string& path) const = 0;
    virtual bool directory_exists(const std::string& path) const = 0;
    virtual FileAttributes get_file_attributes(const std::string& path) const = 0;
    virtual bool set_file_attributes(const std::string& path, const FileAttributes& attributes) = 0;

    // File I/O
    virtual std::unique_ptr<IFileHandle> open_file(const std::string& path, FileOpenMode mode) = 0;
    virtual bool close_file(std::unique_ptr<IFileHandle> handle) = 0;
    virtual bool delete_file(const std::string& path) = 0;
    virtual bool copy_file(const std::string& from, const std::string& to, bool overwrite = false) = 0;
    virtual bool move_file(const std::string& from, const std::string& to, bool overwrite = false) = 0;
    virtual bool create_hard_link(const std::string& target, const std::string& link) = 0;
    virtual bool create_symbolic_link(const std::string& target, const std::string& link) = 0;

    // Directory operations
    virtual bool create_directory(const std::string& path, bool recursive = true) = 0;
    virtual bool delete_directory(const std::string& path, bool recursive = false) = 0;
    virtual bool copy_directory(const std::string& from, const std::string& to, bool overwrite = false) = 0;
    virtual bool move_directory(const std::string& from, const std::string& to, bool overwrite = false) = 0;
    virtual std::vector<DirectoryEntry> list_directory(const std::string& path, bool recursive = false) const = 0;
    virtual std::vector<DirectoryEntry> find_files(const std::string& directory,
                                                  const std::string& pattern,
                                                  bool recursive = false) const = 0;

    // Convenience functions for spherical coordinate data
    virtual bool save_spherical_data(const std::string& path, const void* data, uint64_t size) = 0;
    virtual bool load_spherical_data(const std::string& path, void* data, uint64_t& size) = 0;
    virtual bool save_spherical_scene(const std::string& path, const std::string& scene_data) = 0;
    virtual std::string load_spherical_scene(const std::string& path) = 0;

    // File watching (platform-dependent)
    using FileWatchCallback = std::function<void(const FileWatchEvent& event)>;
    virtual bool start_file_watch(const std::string& path, FileWatchCallback callback) = 0;
    virtual bool stop_file_watch(const std::string& path) = 0;
    virtual bool is_file_watch_supported() const = 0;

    // Asynchronous operations
    virtual bool read_file_async(const std::string& path,
                               std::function<void(const std::vector<uint8_t>&, bool)> callback) = 0;
    virtual bool write_file_async(const std::string& path,
                                const std::vector<uint8_t>& data,
                                std::function<void(bool)> callback) = 0;

    // Memory mapping (for large files)
    virtual bool map_file(const std::string& path, void*& mapped_data, uint64_t& size) = 0;
    virtual bool unmap_file(void* mapped_data) = 0;
    virtual bool is_memory_mapping_supported() const = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;

    // Platform-specific operations
    virtual void* get_platform_handle() const = 0;
};

/**
 * @brief File system factory interface
 */
class IFileSystemFactory {
public:
    virtual ~IFileSystemFactory() = default;
    virtual std::unique_ptr<IFileSystem> create_file_system() = 0;
};

/**
 * @brief Utility functions for file system operations
 */
class FileSystemUtils {
public:
    static std::string read_text_file(IFileSystem& fs, const std::string& path);
    static bool write_text_file(IFileSystem& fs, const std::string& path, const std::string& content);
    static std::vector<uint8_t> read_binary_file(IFileSystem& fs, const std::string& path);
    static bool write_binary_file(IFileSystem& fs, const std::string& path, const std::vector<uint8_t>& data);
    static bool ensure_directory_exists(IFileSystem& fs, const std::string& path);
    static uint64_t calculate_directory_size(IFileSystem& fs, const std::string& path);
    static bool is_path_safe(const std::string& path);
};

} // namespace infrastructure
} // namespace hsml
