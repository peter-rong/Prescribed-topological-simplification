#include "cubical_mode.h"
#include "tet_mode.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

void printUsage(const char* prog) {
    std::cerr
        << "Prescribed topological simplification — cubical | tet\n\n"
        << "Cubical:\n"
        << "  " << prog << " cubical <input.mrc> <out1> <out2> <adjustment> <dtype>\n"
        << "      Optional: [-a|--ascii] [--cpp_program <path>] (ignored) [--topK <n>]\n"
        << "               [--cavitySkip <n>] [--handleSkip <n>] [--componentSkip <n>]\n\n"
        << "Tet:\n"
        << "  " << prog << " tet <mesh.msh> <alpha_file> <out1> <out2> <adjustment> <dtype>\n"
        << "      Optional: [-a|--ascii] [--cpp_program <path>] (ignored) [--topK <n>]\n"
        << "               [--tet_labels <path>]\n"
        << "               [--cavitySkip <n>] [--handleSkip <n>] [--componentSkip <n>]\n";
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string mode = argv[1];
    if (mode == "-h" || mode == "--help" || mode == "/?") {
        printUsage(argv[0]);
        return 0;
    }

    std::vector<char*> shifted;
    shifted.reserve(static_cast<size_t>(argc));
    shifted.push_back(argv[0]);
    for (int i = 2; i < argc; ++i) {
        shifted.push_back(argv[i]);
    }
    const int subArgc = static_cast<int>(shifted.size());

    if (mode == "cubical" || mode == "cubic") {
        if (subArgc < 6) {
            std::cerr << "Cubical mode needs: <input.mrc> <output_file> <output_file2> <adjustment> <dtype>\n\n";
            printUsage(argv[0]);
            return 1;
        }
        return prescribed_topo::cubical::runCubicalMode(subArgc, shifted.data());
    }

    if (mode == "tet") {
        if (subArgc < 7) {
            std::cerr << "Tet mode needs: <mesh.msh> <alpha_file> <output_file> <output_file2> <adjustment> <dtype>\n\n";
            printUsage(argv[0]);
            return 1;
        }
        return prescribed_topo::tet::runTetMode(subArgc, shifted.data());
    }

    std::cerr << "Unknown mode \"" << mode << "\". Expected cubical or tet.\n\n";
    printUsage(argv[0]);
    return 1;
}
