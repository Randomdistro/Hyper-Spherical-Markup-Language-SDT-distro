#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

namespace hsml::cli {

struct CLIOptions {
    std::string command;
    std::unordered_map<std::string, std::string> flags;
    std::vector<std::string> args;
};

class CLI {
public:
    static CLIOptions parse_args(int argc, char** argv) {
        CLIOptions options;

        if (argc < 2) {
            options.command = "run";
            return options;
        }

        options.command = argv[1];

        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg.substr(0, 2) == "--") {
                // Long flag
                size_t eq_pos = arg.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = arg.substr(2, eq_pos - 2);
                    std::string value = arg.substr(eq_pos + 1);
                    options.flags[key] = value;
                } else {
                    options.flags[arg.substr(2)] = "true";
                }
            } else if (arg[0] == '-') {
                // Short flag
                options.flags[arg.substr(1)] = "true";
            } else {
                // Positional argument
                options.args.push_back(arg);
            }
        }

        return options;
    }

    static void print_help() {
        std::cout << "HSML-SDT CLI - Pure Spherical Runtime\n\n";
        std::cout << "Commands:\n";
        std::cout << "  run [file]         Run an HSML file\n";
        std::cout << "  build              Build HSML project\n";
        std::cout << "  test               Run tests\n";
        std::cout << "  validate-physics   Validate SDT physics integrity\n";
        std::cout << "  help               Show this help message\n\n";
        std::cout << "Flags:\n";
        std::cout << "  --port=N           Set server port (default: 3000)\n";
        std::cout << "  --resonance=N      Set resonance frequency (default: 432Hz)\n";
        std::cout << "  --verbose          Enable verbose output\n";
    }
};

} // namespace hsml::cli
