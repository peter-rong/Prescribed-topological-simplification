#include "boundary_matrix.hpp"
#include "topomincut_runner.hpp"

#include <iostream>
#include <limits>
#include <string>

static void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName
              << " [--matrix <matrix_file>] [--alpha <alpha_file>]"
              << " [--topK <number>] [--core <number>] [--neighborhood <number>]"
              << " [--cavitySkip <number>] [--handleSkip <number>] [--componentSkip <number>]"
              << " [--maxTopologyIterations <number>]"
              << " <matrix_file> <alpha_file>" << std::endl;
    std::cerr << "Either use named arguments (--matrix and --alpha) or positional arguments." << std::endl;
}

int main(int argc, char* argv[]) {
    std::string matrix_file;
    std::string alpha_file;
    topomincut::RunParams params;

    if (argc == 3) {
        matrix_file = argv[1];
        alpha_file = argv[2];
    } else if (argc >= 5) {
        for (int i = 1; i < argc; i += 2) {
            std::string arg = argv[i];
            if (arg == "--matrix") {
                matrix_file = argv[i + 1];
            } else if (arg == "--alpha") {
                alpha_file = argv[i + 1];
            } else if (arg == "--topK") {
                params.topK = std::stoi(argv[i + 1]);
                if (params.topK < 0) {
                    std::cerr << "Error: topK must be non-negative" << std::endl;
                    return 1;
                }
            } else if (arg == "--core") {
                const char* v = argv[i + 1];
                if (std::string(v) == "Infinity" || std::string(v) == "Inf") {
                    params.core = std::numeric_limits<double>::infinity();
                } else if (std::string(v) == "-Infinity" || std::string(v) == "-Inf") {
                    params.core = -std::numeric_limits<double>::infinity();
                } else {
                    params.core = std::stod(v);
                }
            } else if (arg == "--neighborhood") {
                const char* v = argv[i + 1];
                if (std::string(v) == "Infinity" || std::string(v) == "Inf") {
                    params.neighborhood = std::numeric_limits<double>::infinity();
                } else if (std::string(v) == "-Infinity" || std::string(v) == "-Inf") {
                    params.neighborhood = -std::numeric_limits<double>::infinity();
                } else {
                    params.neighborhood = std::stod(v);
                }
            } else if (arg == "--cavitySkip") {
                params.cavitySkip = std::stoi(argv[i + 1]);
            } else if (arg == "--handleSkip") {
                params.handleSkip = std::stoi(argv[i + 1]);
            } else if (arg == "--componentSkip") {
                params.componentSkip = std::stoi(argv[i + 1]);
            } else if (arg == "--maxTopologyIterations") {
                params.maxTopologyIterations = std::stoi(argv[i + 1]);
                if (params.maxTopologyIterations < 1) {
                    std::cerr << "Error: maxTopologyIterations must be >= 1" << std::endl;
                    return 1;
                }
            } else {
                printUsage(argv[0]);
                return 1;
            }
        }
    } else {
        printUsage(argv[0]);
        return 1;
    }

    BoundaryMatrix boundary_matrix;
    if (!boundary_matrix.loadFromFile(matrix_file)) {
        std::cerr << "Failed to load matrix from " << matrix_file << std::endl;
        return 1;
    }
    if (!boundary_matrix.loadAlphas(alpha_file)) {
        std::cerr << "Failed to load alpha values from " << alpha_file << std::endl;
        return 1;
    }

    return topomincut::runFromBoundaryMatrix(boundary_matrix, params, nullptr);
}
