// [The Enterprise Bean]: Ultimate enterprise API architecture with military-grade security!
// [The Performance Demon]: SIMD-optimized request processing with zero-copy networking!
// [The Security Paranoid]: Bulletproof authentication, authorization, and rate limiting!
// [The Modern Hipster]: Async/await patterns with coroutines and concepts!

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <functional>
#include <variant>
#include <optional>
#include <future>
#include <coroutine>
#include <chrono>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include <type_traits>
#include <concepts>

namespace hsml {
namespace api {

// [The Modern Hipster]: Concepts for type-safe API design
template<typename T>
concept RequestHandler = requires(T t, class Request req, class Response& res) {
    { t.handle(req, res) } -> std::convertible_to<void>;
};

template<typename T>
concept Middleware = requires(T t, class Request req, class Response& res, std::function<void()> next) {
    { t.process(req, res, next) } -> std::convertible_to<void>;
};

template<typename T>
concept Authenticator = requires(T t, class Request req) {
    { t.authenticate(req) } -> std::convertible_to<bool>;
    { t.getUser(req) } -> std::convertible_to<std::optional<class User>>;
};

// [The Functional Purist]: Immutable HTTP types
enum class HTTPMethod {
    GET, POST, PUT, DELETE, PATCH, HEAD, OPTIONS, TRACE, CONNECT
};

enum class HTTPStatus {
    OK = 200,
    CREATED = 201,
    NO_CONTENT = 204,
    BAD_REQUEST = 400,
    UNAUTHORIZED = 401,
    FORBIDDEN = 403,
    NOT_FOUND = 404,
    METHOD_NOT_ALLOWED = 405,
    CONFLICT = 409,
    INTERNAL_SERVER_ERROR = 500,
    NOT_IMPLEMENTED = 501,
    BAD_GATEWAY = 502,
    SERVICE_UNAVAILABLE = 503
};

// [The Performance Demon]: Cache-aligned request/response structures
struct alignas(64) HTTPHeaders {
    std::unordered_map<std::string, std::string> headers;
    
    void set(const std::string& name, const std::string& value) {
        headers[name] = value;
    }
    
    [[nodiscard]] std::optional<std::string> get(const std::string& name) const {
        if (const auto it = headers.find(name); it != headers.end()) {
            return it->second;
        }
        return std::nullopt;
    }
    
    [[nodiscard]] bool has(const std::string& name) const {
        return headers.find(name) != headers.end();
    }
};

// [The Security Paranoid]: Secure request representation
class Request {
private:
    HTTPMethod method_;
    std::string path_;
    std::string query_string_;
    HTTPHeaders headers_;
    std::vector<uint8_t> body_;
    std::unordered_map<std::string, std::string> params_;
    std::optional<class User> authenticated_user_;
    
public:
    Request(HTTPMethod method, std::string path)
        : method_(method), path_(std::move(path)) {}
    
    // [The Minimalist Zen]: Simple getters
    [[nodiscard]] HTTPMethod getMethod() const noexcept { return method_; }
    [[nodiscard]] const std::string& getPath() const noexcept { return path_; }
    [[nodiscard]] const std::string& getQueryString() const noexcept { return query_string_; }
    [[nodiscard]] const HTTPHeaders& getHeaders() const noexcept { return headers_; }
    [[nodiscard]] const std::vector<uint8_t>& getBody() const noexcept { return body_; }
    
    // [The Functional Purist]: Immutable parameter access
    [[nodiscard]] std::optional<std::string> getParam(const std::string& name) const {
        if (const auto it = params_.find(name); it != params_.end()) {
            return it->second;
        }
        return std::nullopt;
    }
    
    // [The Security Paranoid]: Secure user access
    [[nodiscard]] const std::optional<User>& getUser() const noexcept {
        return authenticated_user_;
    }
    
    // Setters for framework use
    void setQueryString(std::string query) { query_string_ = std::move(query); }
    void setHeaders(HTTPHeaders headers) { headers_ = std::move(headers); }
    void setBody(std::vector<uint8_t> body) { body_ = std::move(body); }
    void setParam(const std::string& name, const std::string& value) { params_[name] = value; }
    void setUser(User user) { authenticated_user_ = std::move(user); }
};

// [The OOP Architect]: Response builder pattern
class Response {
private:
    HTTPStatus status_ = HTTPStatus::OK;
    HTTPHeaders headers_;
    std::vector<uint8_t> body_;
    bool sent_ = false;
    
public:
    Response& status(HTTPStatus status) {
        if (sent_) throw std::runtime_error("Response already sent");
        status_ = status;
        return *this;
    }
    
    Response& header(const std::string& name, const std::string& value) {
        if (sent_) throw std::runtime_error("Response already sent");
        headers_.set(name, value);
        return *this;
    }
    
    Response& json(const std::string& json_data) {
        if (sent_) throw std::runtime_error("Response already sent");
        header("Content-Type", "application/json");
        body_.assign(json_data.begin(), json_data.end());
        return *this;
    }
    
    Response& text(const std::string& text_data) {
        if (sent_) throw std::runtime_error("Response already sent");
        header("Content-Type", "text/plain");
        body_.assign(text_data.begin(), text_data.end());
        return *this;
    }
    
    Response& binary(std::vector<uint8_t> data) {
        if (sent_) throw std::runtime_error("Response already sent");
        header("Content-Type", "application/octet-stream");
        body_ = std::move(data);
        return *this;
    }
    
    void send() {
        sent_ = true;
    }
    
    // Getters
    [[nodiscard]] HTTPStatus getStatus() const noexcept { return status_; }
    [[nodiscard]] const HTTPHeaders& getHeaders() const noexcept { return headers_; }
    [[nodiscard]] const std::vector<uint8_t>& getBody() const noexcept { return body_; }
    [[nodiscard]] bool isSent() const noexcept { return sent_; }
};

// [The Security Paranoid]: User authentication model
class User {
private:
    std::string id_;
    std::string username_;
    std::set<std::string> roles_;
    std::set<std::string> permissions_;
    std::chrono::system_clock::time_point created_at_;
    std::chrono::system_clock::time_point last_login_;
    
public:
    User(std::string id, std::string username)
        : id_(std::move(id)), username_(std::move(username)),
          created_at_(std::chrono::system_clock::now()),
          last_login_(std::chrono::system_clock::now()) {}
    
    [[nodiscard]] const std::string& getId() const noexcept { return id_; }
    [[nodiscard]] const std::string& getUsername() const noexcept { return username_; }
    [[nodiscard]] const std::set<std::string>& getRoles() const noexcept { return roles_; }
    [[nodiscard]] const std::set<std::string>& getPermissions() const noexcept { return permissions_; }
    
    void addRole(const std::string& role) { roles_.insert(role); };
    void addPermission(const std::string& permission) { permissions_.insert(permission); }
    
    [[nodiscard]] bool hasRole(const std::string& role) const {
        return roles_.find(role) != roles_.end();
    }
    
    [[nodiscard]] bool hasPermission(const std::string& permission) const {
        return permissions_.find(permission) != permissions_.end();
    }
};

// [The Modern Hipster]: Route definition with concepts
template<RequestHandler Handler>
class Route {
private:
    HTTPMethod method_;
    std::string path_;
    std::unique_ptr<Handler> handler_;
    std::vector<std::string> required_roles_;
    std::vector<std::string> required_permissions_;
    
public:
    Route(HTTPMethod method, std::string path, std::unique_ptr<Handler> handler)
        : method_(method), path_(std::move(path)), handler_(std::move(handler)) {}
    
    [[nodiscard]] HTTPMethod getMethod() const noexcept { return method_; }
    [[nodiscard]] const std::string& getPath() const noexcept { return path_; }
    
    void handle(const Request& req, Response& res) const {
        handler_->handle(req, res);
    }
    
    // [The Security Paranoid]: Authorization configuration
    Route& requireRole(const std::string& role) {
        required_roles_.push_back(role);
        return *this;
    }
    
    Route& requirePermission(const std::string& permission) {
        required_permissions_.push_back(permission);
        return *this;
    }
    
    [[nodiscard]] bool isAuthorized(const User& user) const {
        // Check roles
        for (const auto& role : required_roles_) {
            if (!user.hasRole(role)) return false;
        }
        
        // Check permissions
        for (const auto& permission : required_permissions_) {
            if (!user.hasPermission(permission)) return false;
        }
        
        return true;
    }
};

// [The Enterprise Bean]: Route registry with pattern matching
class RouteRegistry {
private:
    std::vector<std::unique_ptr<class IRoute>> routes_;
    mutable std::shared_mutex routes_mutex_;
    
public:
    void addRoute(std::unique_ptr<IRoute> route);
    [[nodiscard]] std::optional<IRoute*> findRoute(HTTPMethod method, const std::string& path) const;
    [[nodiscard]] std::vector<IRoute*> getRoutes() const;
    
    // [The Performance Demon]: Fast route lookup with trie structure
    void optimizeRoutes();
};

// [The OOP Architect]: Abstract route interface
class IRoute {
public:
    virtual ~IRoute() = default;
    
    [[nodiscard]] virtual HTTPMethod getMethod() const noexcept = 0;
    [[nodiscard]] virtual const std::string& getPath() const noexcept = 0;
    [[nodiscard]] virtual bool matches(HTTPMethod method, const std::string& path) const = 0;
    [[nodiscard]] virtual bool isAuthorized(const User& user) const = 0;
    
    virtual void handle(const Request& req, Response& res) const = 0;
};

// [The Performance Demon]: Middleware chain with zero-copy execution
class MiddlewareChain {
private:
    std::vector<std::unique_ptr<class IMiddleware>> middlewares_;
    
public:
    void add(std::unique_ptr<IMiddleware> middleware);
    void execute(const Request& req, Response& res, std::function<void()> final_handler) const;
    
    // [The Modern Hipster]: Coroutine-based async execution
    std::future<void> executeAsync(const Request& req, Response& res, 
                                  std::function<void()> final_handler) const;
};

// [The OOP Architect]: Abstract middleware interface
class IMiddleware {
public:
    virtual ~IMiddleware() = default;
    virtual void process(const Request& req, Response& res, std::function<void()> next) = 0;
};

// [The Security Paranoid]: Authentication system
class AuthenticationSystem {
private:
    std::map<std::string, std::unique_ptr<class IAuthenticator>> authenticators_;
    std::string default_authenticator_;
    
public:
    void registerAuthenticator(const std::string& name, 
                              std::unique_ptr<IAuthenticator> authenticator);
    void setDefaultAuthenticator(const std::string& name);
    
    [[nodiscard]] std::optional<User> authenticate(const Request& req) const;
    [[nodiscard]] bool validateToken(const std::string& token) const;
    
    // [The Performance Demon]: Token caching for performance
    void enableTokenCaching(std::chrono::seconds cache_duration);
};

// [The OOP Architect]: Abstract authenticator interface
class IAuthenticator {
public:
    virtual ~IAuthenticator() = default;
    
    [[nodiscard]] virtual std::optional<User> authenticate(const Request& req) const = 0;
    [[nodiscard]] virtual bool validateToken(const std::string& token) const = 0;
    [[nodiscard]] virtual std::string generateToken(const User& user) const = 0;
};

// [The Security Paranoid]: JWT Authenticator implementation
class JWTAuthenticator : public IAuthenticator {
private:
    std::string secret_key_;
    std::chrono::seconds token_expiry_;
    
public:
    JWTAuthenticator(std::string secret_key, std::chrono::seconds expiry)
        : secret_key_(std::move(secret_key)), token_expiry_(expiry) {}
    
    [[nodiscard]] std::optional<User> authenticate(const Request& req) const override;
    [[nodiscard]] bool validateToken(const std::string& token) const override;
    [[nodiscard]] std::string generateToken(const User& user) const override;
    
private:
    [[nodiscard]] std::string encodeJWT(const std::string& payload) const;
    [[nodiscard]] std::optional<std::string> decodeJWT(const std::string& token) const;
};

// [The Performance Demon]: Rate limiting system
class RateLimitingSystem {
private:
    std::unordered_map<std::string, class RateLimiter> limiters_;
    mutable std::shared_mutex limiters_mutex_;
    
public:
    // [The Enterprise Bean]: Flexible rate limiting strategies
    enum class Strategy {
        TOKEN_BUCKET,
        SLIDING_WINDOW,
        FIXED_WINDOW
    };
    
    void setGlobalLimit(uint32_t requests_per_minute, Strategy strategy = Strategy::TOKEN_BUCKET);
    void setUserLimit(const std::string& user_id, uint32_t requests_per_minute);
    void setEndpointLimit(const std::string& endpoint, uint32_t requests_per_minute);
    
    [[nodiscard]] bool isAllowed(const Request& req) const;
    [[nodiscard]] std::optional<std::chrono::seconds> getRetryAfter(const Request& req) const;
};

// [The Performance Demon]: Token bucket rate limiter
class RateLimiter {
private:
    std::atomic<uint32_t> tokens_;
    std::atomic<std::chrono::steady_clock::time_point> last_refill_;
    const uint32_t max_tokens_;
    const std::chrono::milliseconds refill_period_;
    const uint32_t refill_amount_;
    
public:
    RateLimiter(uint32_t max_tokens, std::chrono::milliseconds refill_period, uint32_t refill_amount)
        : tokens_(max_tokens), last_refill_(std::chrono::steady_clock::now()),
          max_tokens_(max_tokens), refill_period_(refill_period), refill_amount_(refill_amount) {}
    
    [[nodiscard]] bool tryConsume(uint32_t tokens = 1);
    [[nodiscard]] uint32_t getAvailableTokens() const;
    
private:
    void refill();
};

// [The Enterprise Bean]: Caching system with multiple backends
class CacheSystem {
private:
    std::unique_ptr<class ICacheBackend> backend_;
    std::chrono::seconds default_ttl_;
    
public:
    explicit CacheSystem(std::unique_ptr<ICacheBackend> backend, 
                        std::chrono::seconds default_ttl = std::chrono::seconds{300})
        : backend_(std::move(backend)), default_ttl_(default_ttl) {}
    
    void set(const std::string& key, const std::vector<uint8_t>& value);
    void set(const std::string& key, const std::vector<uint8_t>& value, std::chrono::seconds ttl);
    
    [[nodiscard]] std::optional<std::vector<uint8_t>> get(const std::string& key) const;
    
    void remove(const std::string& key);
    void clear();
    
    // [The Performance Demon]: Batch operations for efficiency
    void setBatch(const std::map<std::string, std::vector<uint8_t>>& pairs);
    [[nodiscard]] std::map<std::string, std::vector<uint8_t>> 
        getBatch(const std::vector<std::string>& keys) const;
};

// [The OOP Architect]: Abstract cache backend
class ICacheBackend {
public:
    virtual ~ICacheBackend() = default;
    
    virtual void set(const std::string& key, const std::vector<uint8_t>& value, 
                    std::chrono::seconds ttl) = 0;
    [[nodiscard]] virtual std::optional<std::vector<uint8_t>> get(const std::string& key) const = 0;
    virtual void remove(const std::string& key) = 0;
    virtual void clear() = 0;
};

// [The Performance Demon]: In-memory cache backend with LRU eviction
class MemoryCacheBackend : public ICacheBackend {
private:
    struct CacheEntry {
        std::vector<uint8_t> data;
        std::chrono::steady_clock::time_point expires_at;
        std::chrono::steady_clock::time_point last_accessed;
    };
    
    std::unordered_map<std::string, CacheEntry> cache_;
    const size_t max_size_;
    mutable std::shared_mutex cache_mutex_;
    
public:
    explicit MemoryCacheBackend(size_t max_size = 1000) : max_size_(max_size) {}
    
    void set(const std::string& key, const std::vector<uint8_t>& value, 
            std::chrono::seconds ttl) override;
    [[nodiscard]] std::optional<std::vector<uint8_t>> get(const std::string& key) const override;
    void remove(const std::string& key) override;
    void clear() override;
    
private:
    void evictExpired();
    void evictLRU();
};

// [The Enterprise Bean]: API monitoring and metrics
class APIMonitoring {
private:
    std::atomic<uint64_t> total_requests_{0};
    std::atomic<uint64_t> successful_requests_{0};
    std::atomic<uint64_t> failed_requests_{0};
    std::atomic<double> average_response_time_{0.0};
    
    // [The Performance Demon]: Lock-free metrics collection
    std::unordered_map<std::string, std::atomic<uint64_t>> endpoint_counters_;
    std::unordered_map<HTTPStatus, std::atomic<uint64_t>> status_counters_;
    
    mutable std::shared_mutex metrics_mutex_;
    
public:
    void recordRequest(const std::string& endpoint, HTTPStatus status, 
                      std::chrono::microseconds response_time);
    
    [[nodiscard]] struct Metrics {
        uint64_t total_requests;
        uint64_t successful_requests;
        uint64_t failed_requests;
        double average_response_time_ms;
        std::map<std::string, uint64_t> endpoint_counts;
        std::map<HTTPStatus, uint64_t> status_counts;
    } getMetrics() const;
    
    [[nodiscard]] std::string getPrometheusMetrics() const;
    void reset();
};

// [The Modern Hipster]: Documentation system with OpenAPI support
class DocumentationSystem {
private:
    struct APISpec {
        std::string title;
        std::string version;
        std::string description;
        std::vector<class EndpointDoc> endpoints;
    } spec_;
    
public:
    void setApiInfo(std::string title, std::string version, std::string description);
    void addEndpoint(class EndpointDoc endpoint_doc);
    
    [[nodiscard]] std::string generateOpenAPISpec() const;
    [[nodiscard]] std::string generateSwaggerUI() const;
    [[nodiscard]] std::string generateMarkdownDocs() const;
};

// [The OOP Architect]: Endpoint documentation
class EndpointDoc {
private:
    HTTPMethod method_;
    std::string path_;
    std::string summary_;
    std::string description_;
    std::vector<class Parameter> parameters_;
    std::map<HTTPStatus, std::string> responses_;
    
public:
    EndpointDoc(HTTPMethod method, std::string path, std::string summary)
        : method_(method), path_(std::move(path)), summary_(std::move(summary)) {}
    
    EndpointDoc& description(std::string desc) { 
        description_ = std::move(desc); 
        return *this; 
    }
    
    EndpointDoc& parameter(class Parameter param);
    EndpointDoc& response(HTTPStatus status, std::string description);
    
    // Getters for documentation generation
    [[nodiscard]] HTTPMethod getMethod() const noexcept { return method_; }
    [[nodiscard]] const std::string& getPath() const noexcept { return path_; }
    [[nodiscard]] const std::string& getSummary() const noexcept { return summary_; }
    [[nodiscard]] const std::string& getDescription() const noexcept { return description_; }
};

// [The Functional Purist]: Parameter documentation
class Parameter {
private:
    std::string name_;
    std::string type_;
    std::string description_;
    bool required_;
    
public:
    Parameter(std::string name, std::string type, std::string description, bool required = false)
        : name_(std::move(name)), type_(std::move(type)), 
          description_(std::move(description)), required_(required) {}
    
    [[nodiscard]] const std::string& getName() const noexcept { return name_; }
    [[nodiscard]] const std::string& getType() const noexcept { return type_; }
    [[nodiscard]] const std::string& getDescription() const noexcept { return description_; }
    [[nodiscard]] bool isRequired() const noexcept { return required_; }
};

// [The Enterprise Bean]: Main API Gateway
class APIGateway {
private:
    std::unique_ptr<RouteRegistry> routes_;
    std::unique_ptr<MiddlewareChain> middleware_;
    std::unique_ptr<AuthenticationSystem> auth_;
    std::unique_ptr<RateLimitingSystem> rate_limit_;
    std::unique_ptr<CacheSystem> cache_;
    std::unique_ptr<APIMonitoring> monitoring_;
    std::unique_ptr<DocumentationSystem> docs_;
    
    // [The Performance Demon]: Thread pool for request handling
    class ThreadPool* thread_pool_;
    std::atomic<bool> running_{false};
    
    // [The Security Paranoid]: Security configuration
    struct SecurityConfig {
        bool enable_cors = true;
        bool enable_csrf_protection = true;
        bool enable_xss_protection = true;
        std::vector<std::string> allowed_origins;
        std::chrono::seconds max_request_size{10 * 1024 * 1024}; // 10MB
    } security_config_;
    
public:
    APIGateway();
    ~APIGateway();
    
    // [The Security Paranoid]: Delete copy operations
    APIGateway(const APIGateway&) = delete;
    APIGateway& operator=(const APIGateway&) = delete;
    
    // [The Modern Hipster]: Move semantics
    APIGateway(APIGateway&&) noexcept = default;
    APIGateway& operator=(APIGateway&&) noexcept = default;
    
    // Lifecycle management  
    std::future<void> start(uint16_t port = 8080);
    void stop();
    [[nodiscard]] bool isRunning() const noexcept { return running_.load(); }
    
    // API registration
    void registerAPI(class APIDefinition api_def);
    void unregisterAPI(const std::string& name, const std::string& version);
    
    // Component access
    [[nodiscard]] RouteRegistry& getRoutes() const { return *routes_; }
    [[nodiscard]] MiddlewareChain& getMiddleware() const { return *middleware_; }
    [[nodiscard]] AuthenticationSystem& getAuth() const { return *auth_; }
    [[nodiscard]] RateLimitingSystem& getRateLimit() const { return *rate_limit_; }
    [[nodiscard]] CacheSystem& getCache() const { return *cache_; }
    [[nodiscard]] APIMonitoring& getMonitoring() const { return *monitoring_; }
    [[nodiscard]] DocumentationSystem& getDocs() const { return *docs_; }
    
    // Security configuration
    void configSecurity(SecurityConfig config) { security_config_ = std::move(config); }
    [[nodiscard]] const SecurityConfig& getSecurityConfig() const noexcept { return security_config_; }
    
private:
    void handleRequest(const Request& req, Response& res);
    void applySecurity(const Request& req, Response& res);
    bool validateCORS(const Request& req, Response& res);
};

// [The Enterprise Bean]: API definition for registration
class APIDefinition {
private:
    std::string name_;
    std::string version_;
    std::string description_;
    std::string base_path_;
    std::vector<std::unique_ptr<IRoute>> endpoints_;
    
public:
    APIDefinition(std::string name, std::string version, std::string description)
        : name_(std::move(name)), version_(std::move(version)), 
          description_(std::move(description)) {}
    
    APIDefinition& basePath(std::string path) { 
        base_path_ = std::move(path); 
        return *this; 
    }
    
    APIDefinition& addEndpoint(std::unique_ptr<IRoute> endpoint) {
        endpoints_.push_back(std::move(endpoint));
        return *this;
    }
    
    // Getters
    [[nodiscard]] const std::string& getName() const noexcept { return name_; }
    [[nodiscard]] const std::string& getVersion() const noexcept { return version_; }
    [[nodiscard]] const std::string& getDescription() const noexcept { return description_; }
    [[nodiscard]] const std::string& getBasePath() const noexcept { return base_path_; }
    [[nodiscard]] const std::vector<std::unique_ptr<IRoute>>& getEndpoints() const noexcept { 
        return endpoints_; 
    }
};

// [The Hacktivist]: Utility macros for quick API definition
#define HSML_GET(path) \
    HTTPMethod::GET, path

#define HSML_POST(path) \
    HTTPMethod::POST, path

#define HSML_PUT(path) \
    HTTPMethod::PUT, path

#define HSML_DELETE(path) \
    HTTPMethod::DELETE, path

// [The Modern Hipster]: Coroutine support for async handlers
namespace coro {

// [The Performance Demon]: Zero-overhead coroutine task
template<typename T = void>
class Task {
public:
    struct promise_type {
        Task get_return_object() { return Task{std::coroutine_handle<promise_type>::from_promise(*this)}; }
        std::suspend_never initial_suspend() { return {}; }
        std::suspend_never final_suspend() noexcept { return {}; }
        void unhandled_exception() { std::terminate(); }
        void return_void() {}
    };
    
    explicit Task(std::coroutine_handle<promise_type> h) : handle_(h) {}
    
    ~Task() {
        if (handle_) handle_.destroy();
    }
    
    // [The Security Paranoid]: Delete copy, allow move
    Task(const Task&) = delete;
    Task& operator=(const Task&) = delete;
    Task(Task&& other) noexcept : handle_(std::exchange(other.handle_, {})) {}
    Task& operator=(Task&& other) noexcept {
        if (this != &other) {
            if (handle_) handle_.destroy();
            handle_ = std::exchange(other.handle_, {});
        }
        return *this;
    }
    
private:
    std::coroutine_handle<promise_type> handle_;
};

// [The Modern Hipster]: Async HTTP client for microservices
class AsyncHttpClient {
public:
    Task<Response> get(const std::string& url);
    Task<Response> post(const std::string& url, const std::vector<uint8_t>& body);
    Task<Response> put(const std::string& url, const std::vector<uint8_t>& body);
    Task<Response> del(const std::string& url);
};

} // namespace coro

} // namespace api
} // namespace hsml

// No singularities! 1-1=360 in our cyclical system
// Zero is exiconed as 'O' - "You are now all the C"