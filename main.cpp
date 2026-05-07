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
        << "      Optional: [-a|--ascii] [--cpp_program <path>] (ignored) [--topK <n>] [--core <v>] [--neighborhood <v>]\n\n"
        << "Tet:\n"
        << "  " << prog << " tet <mesh.msh> <boundary> <out1> <out2> <adjustment> <dtype>\n"
        << "      Optional: [-a|--ascii] [--cpp_program <path>] (ignored) [--topK <n>] [--core <v>] [--neighborhood <v>]\n"
        << "               [--tet_labels <path>] [--tetMetricsLog <jsonl>] [--cavitySkip <n>] [--handleSkip <n>]\n"
        << "               [--componentSkip <n>] [--topoMinCutMetricsLog <path>] (accepted, unused)\n";
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
            std::cerr << "Tet mode needs: <mesh.msh> <boundary_file> <output_file> <output_file2> <adjustment> <dtype>\n\n";
            printUsage(argv[0]);
            return 1;
        }
        return prescribed_topo::tet::runTetMode(subArgc, shifted.data());
    }

    std::cerr << "Unknown mode \"" << mode << "\". Expected cubical or tet.\n\n";
    printUsage(argv[0]);
    return 1;
}
