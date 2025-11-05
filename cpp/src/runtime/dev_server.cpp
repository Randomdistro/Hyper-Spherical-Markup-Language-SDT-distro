#include <string>
#include <iostream>

namespace hsml::server {

class DevServer {
private:
    int port_;
    bool running_ = false;

public:
    explicit DevServer(int port = 3000) : port_(port) {}

    bool start() {
        std::cout << "Starting HSML-SDT development server on port " << port_ << std::endl;
        running_ = true;
        // Server implementation would go here
        return true;
    }

    void stop() {
        if (running_) {
            std::cout << "Stopping development server" << std::endl;
            running_ = false;
        }
    }

    void handle_request(const std::string& path) {
        // Request handling would go here
        std::cout << "Request: " << path << std::endl;
    }

    bool is_running() const {
        return running_;
    }

    int port() const {
        return port_;
    }
};

} // namespace hsml::server
