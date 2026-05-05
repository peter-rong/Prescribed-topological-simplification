#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <sstream>
#include <iomanip>
#include <limits>
#include <chrono>
#include <cmath>
#include <random>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <tuple>
#include <queue>
#include <numeric>
// Parallel std::sort needs <execution>. On some Linux setups, link with -ltbb or define TOPOTET_NO_PAR_SORT.
#if __has_include(<execution>) && !defined(TOPOTET_NO_PAR_SORT)
#include <execution>
#define TOPOTET_HAS_STD_EXECUTION 1
#endif
#ifdef _WIN32
#include <io.h>
#endif
#include "msh_reader.h"

#include <Eigen/Sparse>
#include "topomincut_runner.hpp"
// Removed marching_cubes.h - not needed for tet mesh
// Removed CGAL includes - alpha values are now read directly from input file

namespace prescribed_topo {
namespace tet {

// TriangleMesh struct for mesh manifoldness checking (separate from Triangle class for cell complex)
struct TriangleMesh {
    int v1, v2, v3;
    TriangleMesh(int v1 = 0, int v2 = 0, int v3 = 0) : v1(v1), v2(v2), v3(v3) {}
};

// Simple manifoldness checker for triangle meshes (edge-based)
struct UndirectedEdgeKey {
    int a;
    int b;
    UndirectedEdgeKey(int u, int v) {
        if (u < v) { a = u; b = v; } else { a = v; b = u; }
    }
    bool operator==(const UndirectedEdgeKey& other) const {
        return a == other.a && b == other.b;
    }
};
struct UndirectedEdgeKeyHash {
    std::size_t operator()(const UndirectedEdgeKey& k) const {
        return (static_cast<std::size_t>(k.a) * 73856093u) ^ (static_cast<std::size_t>(k.b) * 19349663u);
    }
};
static void reportMeshManifoldness(const std::vector<Point3D>& vertices,
                                   const std::vector<TriangleMesh>& triangles,
                                   const std::string& label) {
    std::unordered_map<UndirectedEdgeKey, int, UndirectedEdgeKeyHash> edgeUseCount;
    edgeUseCount.reserve(triangles.size() * 3);
    auto addEdge = [&](int v1, int v2) {
        UndirectedEdgeKey key(v1, v2);
        edgeUseCount[key] += 1;
    };
    for (const auto& t : triangles) {
        // Triangles use 1-based vertex indices in marching_cubes.cpp
        addEdge(t.v1, t.v2);
        addEdge(t.v2, t.v3);
        addEdge(t.v3, t.v1);
    }
    long long totalEdges = 0;
    long long manifoldEdges = 0;
    long long boundaryEdges = 0;     // count == 1
    long long nonManifoldEdges = 0;  // count != 1 && count != 2
    std::vector<std::tuple<int,int,int>> samplesBoundary;
    std::vector<std::tuple<int,int,int>> samplesNonManifold;
    samplesBoundary.reserve(10);
    samplesNonManifold.reserve(10);
    for (const auto& kv : edgeUseCount) {
        totalEdges++;
        int c = kv.second;
        if (c == 2) {
            manifoldEdges++;
        } else if (c == 1) {
            boundaryEdges++;
            if (samplesBoundary.size() < 10) samplesBoundary.emplace_back(kv.first.a, kv.first.b, c);
        } else {
            nonManifoldEdges++;
            if (samplesNonManifold.size() < 10) samplesNonManifold.emplace_back(kv.first.a, kv.first.b, c);
        }
    }
    std::cout << "[Manifoldness] " << label << "\n";
    std::cout << "  Vertices: " << vertices.size() << ", Triangles: " << triangles.size() << "\n";
    std::cout << "  Unique undirected edges: " << totalEdges << "\n";
    std::cout << "  Manifold edges (count==2): " << manifoldEdges << "\n";
    std::cout << "  Boundary edges (count==1): " << boundaryEdges << "\n";
    std::cout << "  Non-manifold edges (count>2): " << nonManifoldEdges << "\n";
    if (!samplesBoundary.empty()) {
        std::cout << "  Sample boundary edges (v1,v2,count): ";
        for (size_t i = 0; i < samplesBoundary.size(); ++i) {
            auto [a,b,c] = samplesBoundary[i];
            std::cout << "(" << a << "," << b << "," << c << ")";
            //also print the vertex coords
            std::cout << " (" << vertices[a].x << "," << vertices[a].y << "," << vertices[a].z << ")-(" << vertices[b].x << "," << vertices[b].y << "," << vertices[b].z << ")";
            if (i + 1 < samplesBoundary.size()) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    if (!samplesNonManifold.empty()) {
        std::cout << "  Sample non-manifold edges (v1,v2,count): ";
        for (size_t i = 0; i < samplesNonManifold.size(); ++i) {
            auto [a,b,c] = samplesNonManifold[i];
            std::cout << "(" << a << "," << b << "," << c << ")";
            if (i + 1 < samplesNonManifold.size()) std::cout << ", ";
        }
        std::cout << std::endl;
    }
}

// Forward declarations
class Simplex;
class Vertex;
class Edge;
class Triangle;
class Tet;


// Base class for all simplices (memory optimized)
class Simplex {
public:
    int id;
    double val;
    int dimension;
    double x,y,z;
    bool changed = false;
    std::vector<Simplex*> faces;
    
    virtual ~Simplex() = default;
    
    bool operator<(const Simplex& other) const {
        if (val == other.val) {
            if (dimension == other.dimension) {
                return id < other.id;
            } else {
                return dimension < other.dimension;
            }
        } else {
            return val < other.val;
        }
    }
    
    // Return by reference — returning a copy of `faces` was O(faces) per call and dominated
    // childrenByColumn construction (millions of allocations on large meshes).
    virtual const std::vector<Simplex*>& getFaces() const {
        return faces;
    }
    
    // Memory optimization: shrink faces vector after construction
    virtual void optimizeMemory() {
        faces.shrink_to_fit();
    }
};

class Vertex : public Simplex {
public:
    std::vector<Edge*> incident_edges;
    
    Vertex(int id, double val, double x, double y, double z) {
        this->id = id;
        this->val = val;
        this->dimension = 0;
        this->x = x;
        this->y = y;
        this->z = z;
        this->changed = false;
        // Pre-allocate incident_edges to avoid reallocations
        incident_edges.reserve(12); // Max 12 edges per vertex in 3D grid
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
        incident_edges.shrink_to_fit();
    }
};

class Edge : public Simplex {
public:
    std::vector<Triangle*> incident_triangles;
    std::array<std::array<double, 3>, 2> verts{};  // Fixed-size endpoints [2][3]
    
    Edge(int id, double val, const std::array<std::array<double, 3>, 2>& verts, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 1;
        this->faces = children;
        this->verts = verts;
        this->changed = false;
        // Pre-allocate incident_triangles
        incident_triangles.reserve(6); // Max triangles per edge in tet mesh (can be shared by multiple tets)
    }

    Edge(int id, double val, std::array<std::array<double, 3>, 2>&& verts, std::vector<Simplex*>&& children) {
        this->id = id;
        this->val = val;
        this->dimension = 1;
        this->faces = std::move(children);
        this->verts = std::move(verts);
        this->changed = false;
        incident_triangles.reserve(6);
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
        incident_triangles.shrink_to_fit();
    }
};

class Triangle : public Simplex {
public:
    std::vector<Tet*> incident_tets;
    std::array<std::array<double, 3>, 3> verts{};  // Fixed-size vertices [3][3]
    
    Triangle(int id, double val, const std::array<std::array<double, 3>, 3>& verts, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 2;
        this->faces = children;
        this->verts = verts;
        // Pre-allocate incident_tets
        incident_tets.reserve(2); // Max 2 tets per triangle (boundary triangles have 1)
    }

    Triangle(int id, double val, std::array<std::array<double, 3>, 3>&& verts, std::vector<Simplex*>&& children) {
        this->id = id;
        this->val = val;
        this->dimension = 2;
        this->faces = std::move(children);
        this->verts = std::move(verts);
        incident_tets.reserve(2);
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
        incident_tets.shrink_to_fit();
    }
};

class Tet : public Simplex {
public:
    std::array<std::array<double, 3>, 4> verts{};  // Fixed-size vertices [4][3]
    
    Tet(int id, double val, const std::array<std::array<double, 3>, 4>& verts, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 3;
        this->faces = children;
        this->verts = verts;
        this->changed = false;
    }

    Tet(int id, double val, std::array<std::array<double, 3>, 4>&& verts, std::vector<Simplex*>&& children) {
        this->id = id;
        this->val = val;
        this->dimension = 3;
        this->faces = std::move(children);
        this->verts = std::move(verts);
        this->changed = false;
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
    }
};

namespace {

// Inline sort key (same order as Simplex::operator<) — used in std::sort hot path.
inline bool simplexPtrLess(const Simplex* a, const Simplex* b) {
    const double av = a->val;
    const double bv = b->val;
    if (av != bv) {
        return av < bv;
    }
    const int ad = a->dimension;
    const int bd = b->dimension;
    if (ad != bd) {
        return ad < bd;
    }
    return a->id < b->id;
}

} // namespace

// Filter result for keeping only in-range ∪ interface simplices (memory optimized)
struct FilterResult {
    std::vector<int> keptDims;                         // kept simplex dimensions in new order
    std::vector<double> keptVals;                      // kept simplex values in new order
    std::vector<std::vector<int>> keptChildrenByCol;   // child rows per kept column (0-based), new indices
    std::vector<int> oldToNew;                         // size N, -1 if dropped
    std::vector<int> newToOld;                         // size K
    
    // Memory optimization: shrink vectors after construction
    void optimizeMemory() {
        keptDims.shrink_to_fit();
        keptVals.shrink_to_fit();
        for (auto& col : keptChildrenByCol) {
            col.shrink_to_fit();
        }
        keptChildrenByCol.shrink_to_fit();
        oldToNew.shrink_to_fit();
        newToOld.shrink_to_fit();
    }
};

// Filter to keep inRange (alpha1 <= val <= alpha2) and interface (neighbor of inRange) simplices
inline FilterResult filterInterfaceAndInRange(
    const std::vector<int>& simplexDims,
    const std::vector<double>& simplexVals,
    const std::vector<std::vector<int>>& childrenByColumn,
    double alpha1,
    double alpha2
    )
{
    const int n = static_cast<int>(simplexDims.size());
    std::vector<bool> inRange(n, false); // Use bool instead of char for memory efficiency
    for (int i = 0; i < n; ++i) {
        const double a = simplexVals[static_cast<size_t>(i)];
        (void)alpha1;
        (void)alpha2;
        (void)a;
        inRange[i] = true;
    }

    // Fast path: every cell inRange → keep all, identity map (skip O(E) parents transpose).
    bool allInRange = true;
    for (int i = 0; i < n; ++i) {
        if (!inRange[i]) {
            allInRange = false;
            break;
        }
    }
    if (allInRange && n > 0) {
        FilterResult result;
        result.keptDims = simplexDims;
        result.keptVals = simplexVals;
        result.keptChildrenByCol = childrenByColumn;
        result.oldToNew.resize(n);
        result.newToOld.resize(n);
        for (int i = 0; i < n; ++i) {
            result.oldToNew[i] = i;
            result.newToOld[i] = i;
        }
        return result;
    }

    //print how many are inRange for each dimension
    int totalInRangeVertices = 0;
    int totalInRangeEdges = 0;
    int totalInRangeTriangles = 0;
    int totalInRangeTets = 0;
    for (int i = 0; i < n; ++i) {
        if (inRange[i]){
            if (simplexDims[static_cast<size_t>(i)] == 0) totalInRangeVertices++;
            if (simplexDims[static_cast<size_t>(i)] == 1) totalInRangeEdges++;
            if (simplexDims[static_cast<size_t>(i)] == 2) totalInRangeTriangles++;
            if (simplexDims[static_cast<size_t>(i)] == 3) totalInRangeTets++;
        }
    }
    std::cout << "Total inRange vertices: " << totalInRangeVertices << std::endl;
    std::cout << "Total inRange edges: " << totalInRangeEdges << std::endl;
    std::cout << "Total inRange triangles: " << totalInRangeTriangles << std::endl;
    std::cout << "Total inRange tets: " << totalInRangeTets << std::endl;

    //find the first and last index of the inRange array of 1
    int firstIndex = -1;
    int lastIndex = -1;
    for (int i = 0; i < n; ++i) {
        if (inRange[i]) {
            if (firstIndex == -1) firstIndex = i;
            lastIndex = i;
        }
    }
    std::cout << n << std::endl;

    // Build parents (transpose childrenByColumn)
    std::vector<std::vector<int>> parentsByRow(n);
    parentsByRow.assign(n, {});
    
    for (int j = 0; j < n; ++j) {
        for (int child : childrenByColumn[j]) {
            if (child >= 0 && child < n) parentsByRow[child].push_back(j);
        }
    }

    // Interface if any parent or child is inRange
    std::vector<bool> isInterface(n, false); // Use bool instead of char for memory efficiency

    for (int i = firstIndex; i > 0; --i) {
        if (inRange[i]) continue;
        for (int p : parentsByRow[i]) { if (inRange[p]||isInterface[p]) { isInterface[i] = true; break; }}

    }

    //output how many are interface
    int totalInterface = 0;
    int totalInterfaceVertices = 0;
    int totalInterfaceEdges = 0;
    int totalInterfaceTriangles = 0;
    int totalInterfaceTets = 0;

    for (int i = 0; i < isInterface.size(); ++i) {
        if (isInterface[i]){ 
            
            totalInterface++;

            if (simplexDims[static_cast<size_t>(i)] == 0) totalInterfaceVertices++;
            if (simplexDims[static_cast<size_t>(i)] == 1) totalInterfaceEdges++;
            if (simplexDims[static_cast<size_t>(i)] == 2) totalInterfaceTriangles++;
            if (simplexDims[static_cast<size_t>(i)] == 3) totalInterfaceTets++;
        }
    }

    for (int i = lastIndex; i < n; ++i) {
        if (inRange[i]) continue;
        for (int c : childrenByColumn[i]) { if (inRange[c]||isInterface[c]) { isInterface[i] = true; break; } }

    }

    totalInterface = 0;
    totalInterfaceVertices = 0;
    totalInterfaceEdges = 0;
    totalInterfaceTriangles = 0;
    totalInterfaceTets = 0;

    for (int i = 0; i < isInterface.size(); ++i) {
        if (isInterface[i]){ 
            
            totalInterface++;

            if (simplexDims[static_cast<size_t>(i)] == 0) totalInterfaceVertices++;
            if (simplexDims[static_cast<size_t>(i)] == 1) totalInterfaceEdges++;
            if (simplexDims[static_cast<size_t>(i)] == 2) totalInterfaceTriangles++;
            if (simplexDims[static_cast<size_t>(i)] == 3) totalInterfaceTets++;
        }
    }

    // Keep mapping
    std::vector<int> oldToNew(n, -1);
    std::vector<int> newToOld;
    newToOld.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (inRange[i] || isInterface[i]) {
            oldToNew[i] = static_cast<int>(newToOld.size());
            newToOld.push_back(i);
        }
    }
    const int k = static_cast<int>(newToOld.size());

    std::vector<int> keptDims;
    std::vector<double> keptVals;
    keptDims.reserve(static_cast<size_t>(k));
    keptVals.reserve(static_cast<size_t>(k));
    for (int idx : newToOld) {
        keptDims.push_back(simplexDims[static_cast<size_t>(idx)]);
        keptVals.push_back(simplexVals[static_cast<size_t>(idx)]);
    }

    // Kept boundary rows per kept column
    std::vector<std::vector<int>> keptChildrenByCol(k);
    for (int newCol = 0; newCol < k; ++newCol) {
        int oldCol = newToOld[newCol];
        const auto& oldChildren = childrenByColumn[oldCol];
        auto& newChildren = keptChildrenByCol[newCol];
        newChildren.reserve(oldChildren.size());
        for (int oldRow : oldChildren) {
            int newRow = (oldRow >= 0 && oldRow < n) ? oldToNew[oldRow] : -1;
            if (newRow >= 0) newChildren.push_back(newRow);
        }
        // Optional dedup/sort
        /*
        std::sort(newChildren.begin(), newChildren.end());
        newChildren.erase(std::unique(newChildren.begin(), newChildren.end()), newChildren.end());
        */
    }

    FilterResult result{ std::move(keptDims), std::move(keptVals), std::move(keptChildrenByCol), std::move(oldToNew), std::move(newToOld) };
    result.optimizeMemory();
    return result;
}

// Helper function to read cell indices from file
std::vector<int> read_cell_indices(const std::string& filename) {
    std::vector<int> indices;
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        try {
            int index = std::stoi(line);
            indices.push_back(index);
        } catch (const std::exception& e) {
            continue;
        }
    }
    
    return indices;
}

// Helper function to read double values from file
std::vector<double> read_double_values(const std::string& filename) {
    std::vector<double> values;
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        try {
            double value = std::stod(line);
            values.push_back(value);
        } catch (const std::exception& e) {
            continue;
        }
    }

    return values;
}

// Helper function to execute shell command with real-time output
std::string exec_command(const std::string& cmd) {
#ifdef _WIN32
    // Windows implementation using _popen
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(cmd.c_str(), "r"), _pclose);
    
    if (!pipe) {
        throw std::runtime_error("_popen() failed!");
    }
    
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
        std::cout << buffer.data(); // Print output in real-time
        result += buffer.data();
    }
    
    return result;
#else
    // Unix/Linux implementation using popen
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
    
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) != nullptr) {
        std::cout << buffer.data(); // Print output in real-time
        result += buffer.data();
    }
    
    return result;
#endif
}

// Helper function to convert int to bytes
std::vector<uint8_t> getBytesFromInt(int64_t i) {
    std::vector<uint8_t> bytes(8);
    for (int j = 0; j < 8; ++j) {
        bytes[j] = (i >> (8 * j)) & 0xFF;
    }
    return bytes;
}

inline bool sameSign(double a, double b) {
    return (a > 0.0 && b > 0.0) || (a <= 0.0 && b <= 0.0);
}

// Helper structure for edge indexing (undirected)
struct EdgeKey {
    int v1, v2;
    EdgeKey(int a, int b) {
        if (a < b) { v1 = a; v2 = b; }
        else { v1 = b; v2 = a; }
    }
    bool operator==(const EdgeKey& other) const {
        return v1 == other.v1 && v2 == other.v2;
    }
    bool operator<(const EdgeKey& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const {
        return (static_cast<std::size_t>(k.v1) * 73856093u) ^ (static_cast<std::size_t>(k.v2) * 19349663u);
    }
};

// Helper structure for triangle indexing (undirected, sorted vertices)
struct TriangleKey {
    int v1, v2, v3;
    TriangleKey(int a, int b, int c) {
        // Sort vertices
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        v1 = a; v2 = b; v3 = c;
    }
    bool operator==(const TriangleKey& other) const {
        return v1 == other.v1 && v2 == other.v2 && v3 == other.v3;
    }
    bool operator<(const TriangleKey& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        if (v2 != other.v2) return v2 < other.v2;
        return v3 < other.v3;
    }
};

struct TriangleKeyHash {
    std::size_t operator()(const TriangleKey& k) const {
        return (static_cast<std::size_t>(k.v1) * 73856093u) ^ 
               (static_cast<std::size_t>(k.v2) * 19349663u) ^ 
               (static_cast<std::size_t>(k.v3) * 83492791u);
    }
};

int runTetMode(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] 
                  << " tet <mesh_file> <boundary_file> <output_file> <output_file2> <adjustment> <dtype> "
                  << "[-a] [--cpp_program <path>] [--topK <value>] [--core <value>] [--neighborhood <value>] "
                  << "[--tet_labels <path>] "
                  << "[--tetMetricsLog <jsonl_file>] "
                  << "[--cavitySkip <number>] [--handleSkip <number>] [--componentSkip <number>] "
                  << "[--topoMinCutMetricsLog <jsonl_file>]" << std::endl;
        return 1;
    }
    
    std::string mesh_file = argv[1];
    std::string boundary_file = argv[2];
    std::string output_file = "output/" + std::string(argv[3]);
    std::string output_file2 = "output/" + std::string(argv[4]);
    std::string adjustment_str = argv[5];
    std::string dtype = argv[6];
    
    bool ascii_mode = false;
    std::string cpp_program = "C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\build\\release\\TopoMinCut";
    int topK = 0;
    std::string tet_labels_file = "";  // Optional tet labels file

    // JSONL logging (this driver) and metrics/logging for TopoMinCut.
    std::string tetMetricsLogFile = "";
    std::string topoMinCutMetricsLogFile = "";
    int cavitySkip = 0;
    int handleSkip = 0;
    int componentSkip = 0;

    // core/neighborhood use floating infinities so callers can check via std::isinf
    double core = std::numeric_limits<double>::infinity();
    double neighborhood = -std::numeric_limits<double>::infinity();
    bool create_3d_array = true;
    
    // Parse optional arguments
    for (int i = 7; i < argc; ++i) {
        if (std::string(argv[i]) == "-a" || std::string(argv[i]) == "--ascii") {
            ascii_mode = true;
        } else if (std::string(argv[i]) == "--cpp_program" && i + 1 < argc) {
            cpp_program = argv[++i];
        } else if (std::string(argv[i]) == "--topK" && i + 1 < argc) {
            topK = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--core" && i + 1 < argc) {
            core = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == "--neighborhood" && i + 1 < argc) {
            neighborhood = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == "--tet_labels" && i + 1 < argc) {
            tet_labels_file = argv[++i];
        } else if (std::string(argv[i]) == "--tetMetricsLog" && i + 1 < argc) {
            tetMetricsLogFile = argv[++i];
        } else if (std::string(argv[i]) == "--cavitySkip" && i + 1 < argc) {
            cavitySkip = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--handleSkip" && i + 1 < argc) {
            handleSkip = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--componentSkip" && i + 1 < argc) {
            componentSkip = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--topoMinCutMetricsLog" && i + 1 < argc) {
            topoMinCutMetricsLogFile = argv[++i];
        }
    }
    (void)topoMinCutMetricsLogFile; // TopoMinCut runner no longer writes metrics files

    // Timings to write into JSONL (when enabled).
    double simplex2TimeMs = -1.0;
    double isosurfaceTimeMs = -1.0;

    // Append a single JSON object per run (JSON Lines format).
    auto appendJsonlMetrics = [&](double simplex2Ms, double isosurfaceMs) {
        if (tetMetricsLogFile.empty()) return;

        std::filesystem::path p(tetMetricsLogFile);
        if (!p.parent_path().empty()) {
            std::filesystem::create_directories(p.parent_path());
        }

        std::ofstream out(tetMetricsLogFile, std::ios::out | std::ios::app);
        if (!out) {
            std::cerr << "Warning: Failed to open tet metrics log for appending: " << tetMetricsLogFile << std::endl;
            return;
        }

        auto jsonNumberOrNull = [&](double v) -> std::string {
            if (std::isfinite(v)) {
                std::ostringstream ss;
                ss << std::setprecision(17) << v;
                return ss.str();
            }
            return "null";
        };

        out << "{"
            << "\"simplex2_ms\":" << jsonNumberOrNull(simplex2Ms) << ","
            << "\"isosurface_ms\":" << jsonNumberOrNull(isosurfaceMs)
            << "}\n";
    };

    double adjustment = std::stod(adjustment_str);

    (void)cpp_program; // TopoMinCut runs in-process; flag kept for backward-compatible CLIs.
    
    core = -(core - adjustment);
    neighborhood = -(neighborhood - adjustment);

    std::cout << "Adjusted Core: " << core << std::endl;
    std::cout << "Adjusted Neighborhood: " << neighborhood << std::endl;
    
    // Read mesh file (.msh)
    std::cout << "Reading mesh file: " << mesh_file << std::endl;
    
    std::vector<Point3D> vertices;
    std::vector<TetElement> tetrahedra;
    
    try {
        auto mesh_data = MSHReader::readMSHFile(mesh_file);
        vertices = mesh_data.first;
        tetrahedra = mesh_data.second;
    } catch (const std::exception& e) {
        std::cerr << "Error reading mesh file: " << e.what() << std::endl;
        return 1;
    }
    
    // Read vertex alpha values
    std::cout << "Reading vertex alpha values file: " << boundary_file << std::endl;
    std::vector<double> vertexAlphas;
    try {

        //adjustment is applied to the alpha values in the file
        vertexAlphas = MSHReader::readVertexAlphas(boundary_file, adjustment);
    } catch (const std::exception& e) {
        std::cerr << "Error reading vertex alpha values file: " << e.what() << std::endl;
        return 1;
    }
    
    // Verify we have alpha values for all vertices
    if (vertexAlphas.size() != vertices.size()) {
        std::cerr << "Error: Number of alpha values (" << vertexAlphas.size() 
                  << ") does not match number of vertices (" << vertices.size() << ")" << std::endl;
        return 1;
    }
    /*
    // Read sign file
    std::cout << "Reading sign file: " << sign_file << std::endl;
    std::vector<int> signs;
    try {
        signs = MSHReader::readSignFile(sign_file);
        std::cout << "Read " << signs.size() << " sign entries" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error reading sign file: " << e.what() << std::endl;
        return 1;
    }
    
    // Verify we have signs for all tetrahedra
    if (signs.size() != tetrahedra.size()) {
        std::cerr << "Error: Number of signs (" << signs.size() 
                  << ") does not match number of tetrahedra (" << tetrahedra.size() << ")" << std::endl;
        return 1;
    }
    
    // Adjust vertex alpha values based on tet signs
    // For tets with sign = -1: if any vertex has positive alpha, flip it to negative
    // For tets with sign = +1: if any vertex has negative alpha, flip it to positive
    int changedVertices = 0;  // Track which vertices were changed (avoid double counting)
    
    for (size_t tetIdx = 0; tetIdx < tetrahedra.size(); ++tetIdx) {
        int tetSign = signs[tetIdx];
        const auto& tet = tetrahedra[tetIdx];
        
        // Get all 4 vertices of this tet
        int v1 = tet.v1;
        int v2 = tet.v2;
        int v3 = tet.v3;
        int v4 = tet.v4;
        
        if (tetSign == 1) {
            // For sign = -1: flip positive vertex alphas to negative
            if (v1 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v1] > 0.0) {
                vertexAlphas[v1] *= -1.0;
                changedVertices++;
            }
            if (v2 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v2] > 0.0) {
                vertexAlphas[v2] *= -1.0;
                changedVertices++;
            }
            if (v3 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v3] > 0.0) {
                vertexAlphas[v3] *= -1.0;
                changedVertices++;
            }
            if (v4 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v4] > 0.0) {
                vertexAlphas[v4] *= -1.0;
                changedVertices++;
            }
        } else if (tetSign == -1) {
            // For sign = +1: flip negative vertex alphas to positive
            if (v1 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v1] < 0.0) {
                vertexAlphas[v1] *= -1.0;
                changedVertices++;
            }
            if (v2 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v2] < 0.0) {
                vertexAlphas[v2] *= -1.0;
                changedVertices++;
            }
            if (v3 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v3] < 0.0) {
                vertexAlphas[v3] *= -1.0;
                changedVertices++;
            }
            if (v4 < static_cast<int>(vertexAlphas.size()) && vertexAlphas[v4] < 0.0) {
                vertexAlphas[v4] *= -1.0;
                changedVertices++;
            }
        }
        // For sign = 0, do nothing
    }
    
    std::cout << "Adjusted vertex alpha values: " << changedVertices
              << " unique vertices were flipped based on tet signs" << std::endl;
    */
    std::cout << "Mesh has " << vertices.size() << " vertices and " << tetrahedra.size() << " tetrahedra" << std::endl;

    
    // Read tet labels from tet_labels file (if provided)
    std::vector<int> tetLabels;
    if (!tet_labels_file.empty()) {
        std::cout << "Reading tet labels file: " << tet_labels_file << std::endl;
        try {
            tetLabels = MSHReader::readTetLabels(tet_labels_file);
            if (tetLabels.size() != tetrahedra.size()) {
                std::cerr << "Warning: Number of tet labels (" << tetLabels.size() 
                          << ") does not match number of tetrahedra (" << tetrahedra.size() << ")" << std::endl;
                std::cerr << "Will use default labels (all outside)" << std::endl;
                tetLabels.clear();
                tetLabels.resize(tetrahedra.size(), 0);  // Default to 0 (outside)
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not read tet labels file: " << e.what() << std::endl;
            std::cerr << "Will use default labels (all outside)" << std::endl;
            tetLabels.resize(tetrahedra.size(), 0);  // Default to 0 (outside)
        }
    } else {
        std::cout << "No tet labels file provided. Will use default labels (all outside)" << std::endl;
        tetLabels.resize(tetrahedra.size(), 0);  // Default to 0 (outside)
    }

    // Compute XY bounding box size and perturb vertex alpha values.
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    for (const auto& p : vertices) {
        minX = std::min(minX, p.x);
        maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }

    const double dx = maxX - minX;
    const double dy = maxY - minY;
    // Use the 2D diagonal of the XY bounding box as the scale.
    const double boundingBoxSize = std::sqrt(dx * dx + dy * dy);
    
    const double epsilon = 0;
    //const double epsilon = 0.01;
   
    const double perturbationAmplitude = epsilon * boundingBoxSize;

    
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> perturbDist(-perturbationAmplitude, perturbationAmplitude);
    for (auto& alpha : vertexAlphas) {
        alpha += perturbDist(rng);
    }
    
    // Rebuild tet labels from perturbed vertex alpha signs:
    // -1 if all 4 vertices are non-positive, +1 otherwise.
    tetLabels.resize(tetrahedra.size(), 1);
    for (size_t tetIdx = 0; tetIdx < tetrahedra.size(); ++tetIdx) {
        const auto& tet = tetrahedra[tetIdx];
        const bool allNonPositive =
            (vertexAlphas[tet.v1] <= 0.0) &&
            (vertexAlphas[tet.v2] <= 0.0) &&
            (vertexAlphas[tet.v3] <= 0.0) &&
            (vertexAlphas[tet.v4] <= 0.0);
        tetLabels[tetIdx] = allNonPositive ? -1 : 1;
    }
    
    // ============================================================================
    // MARCHING TETRAHEDRA IMPLEMENTATION (before topoMinCut)
    // ============================================================================
    
    // Tet edge map: edges are defined by vertex pairs (1-indexed in lookup table)
    // tetEdges[0] = {1,2} -> vertices {v1, v2}
    // tetEdges[1] = {1,3} -> vertices {v1, v3}
    // tetEdges[2] = {1,4} -> vertices {v1, v4}
    // tetEdges[3] = {2,3} -> vertices {v2, v3}
    // tetEdges[4] = {2,4} -> vertices {v2, v4}
    // tetEdges[5] = {3,4} -> vertices {v3, v4}
    const std::array<std::pair<int, int>, 6> tetEdges = {{
        {0, 1},  // edge 1: v1-v2 (1-indexed: {1,2})
        {0, 2},  // edge 2: v1-v3 (1-indexed: {1,3})
        {0, 3},  // edge 3: v1-v4 (1-indexed: {1,4})
        {1, 2},  // edge 4: v2-v3 (1-indexed: {2,3})
        {1, 3},  // edge 5: v2-v4 (1-indexed: {2,4})
        {2, 3}   // edge 6: v3-v4 (1-indexed: {3,4})
    }};
    
    // Marching-tet lookup table: each entry is {edge_indices, triangle_cycles}
    // Edge indices are 1-indexed in the table, so we subtract 1 when using them
    // Pattern index = binary pattern from vertex signs (0-15)
    struct TetTableEntry {
        std::vector<int> edgeIndices;      // 1-indexed edge indices from table
        std::vector<std::vector<int>> triangles;  // Triangle cycles (edge indices)
    };
    
    const std::array<TetTableEntry, 16> tetTableEinds = {{
        {{}, {}},  // Pattern 0: all negative
        {{3, 5, 6}, {{3, 5, 6}}},  // Pattern 1: only v4 positive
        {{2, 6, 4}, {{2, 6, 4}}},  // Pattern 2: only v3 positive
        {{2, 3, 5, 4}, {{2, 3, 5}, {2, 5, 4}}},  // Pattern 3: v3, v4 positive
        {{1, 4, 5}, {{1, 4, 5}}},  // Pattern 4: only v2 positive
        {{1, 4, 6, 3}, {{1, 6, 3}, {1, 4, 6}}},  // Pattern 5: v2, v4 positive
        {{1, 2, 6, 5}, {{1, 2, 6}, {1, 6, 5}}},  // Pattern 6: v2, v3 positive
        {{1, 2, 3}, {{1, 2, 3}}},  // Pattern 7: v2, v3, v4 positive
        {{1, 3, 2}, {{1, 3, 2}}},  // Pattern 8: only v1 positive
        {{1, 5, 6, 2}, {{1, 6, 2}, {1, 5, 6}}},  // Pattern 9: v1, v4 positive
        {{1, 3, 6, 4}, {{1, 3, 6}, {1, 6, 4}}},  // Pattern 10: v1, v3 positive
        {{1, 5, 4}, {{1, 5, 4}}},  // Pattern 11: v1, v3, v4 positive
        {{2, 4, 5, 3}, {{2, 5, 3}, {2, 4, 5}}},  // Pattern 12: v1, v2 positive
        {{2, 4, 6}, {{2, 4, 6}}},  // Pattern 13: v1, v2, v4 positive
        {{3, 6, 5}, {{3, 6, 5}}},  // Pattern 14: v1, v2, v3 positive
        {{}, {}}   // Pattern 15: all positive
    }};
    
    // Helper: Get sign (0 for <= 0, 1 for > 0)
    auto getSign = [](double val) -> int {
        return (val > 0.0) ? 1 : 0;
    };
    
    // Helper: Compute tet orientation (returns sign of determinant)
    // getTetOri[{a, b, c, d}] := Sign[Cross[b - a, c - a] . (d - a)]
    auto getTetOrientation = [](const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) -> int {
        // Compute b - a, c - a, d - a
        double ba_x = b.x - a.x, ba_y = b.y - a.y, ba_z = b.z - a.z;
        double ca_x = c.x - a.x, ca_y = c.y - a.y, ca_z = c.z - a.z;
        double da_x = d.x - a.x, da_y = d.y - a.y, da_z = d.z - a.z;
        
        // Cross product: (b - a) × (c - a)
        double cross_x = ba_y * ca_z - ba_z * ca_y;
        double cross_y = ba_z * ca_x - ba_x * ca_z;
        double cross_z = ba_x * ca_y - ba_y * ca_x;
        
        // Dot product with (d - a)
        double dot = cross_x * da_x + cross_y * da_y + cross_z * da_z;
        
        return (dot > 0.0) ? 1 : ((dot < 0.0) ? -1 : 0);
    };
    
    // Helper: Compute zero-crossing point on edge using linear interpolation
    auto getLinearRoot = [](const Point3D& v1, const Point3D& v2, double val1, double val2) -> Point3D {
        // Avoid division by zero - if values are same, return midpoint
        if (std::abs(val1 - val2) < 1e-12) {
            return Point3D((v1.x + v2.x) * 0.5, (v1.y + v2.y) * 0.5, (v1.z + v2.z) * 0.5);
        }
        
        // Linear interpolation: t such that val1 + t*(val2 - val1) = 0
        double t = -val1 / (val2 - val1);
        
        // Clamp t away from endpoints to avoid numerical issues with very small triangles
        // Use epsilon 0.01 to push away from 0 and 1
        const double eps = 0.1;
        if (t < eps) t = eps;
        else if (t > 1.0 - eps) t = 1.0 - eps;
        
        return Point3D(
            v1.x + t * (v2.x - v1.x),
            v1.y + t * (v2.y - v1.y),
            v1.z + t * (v2.z - v1.z)
        );
    };
    
    // Helper: Compute triangle area (squared) to check for near-degenerate triangles
    auto getTriangleAreaSquared = [](const Point3D& v1, const Point3D& v2, const Point3D& v3) -> double {
        // Compute two edge vectors
        double e1x = v2.x - v1.x;
        double e1y = v2.y - v1.y;
        double e1z = v2.z - v1.z;
        
        double e2x = v3.x - v1.x;
        double e2y = v3.y - v1.y;
        double e2z = v3.z - v1.z;
        
        // Cross product
        double cross_x = e1y * e2z - e1z * e2y;
        double cross_y = e1z * e2x - e1x * e2z;
        double cross_z = e1x * e2y - e1y * e2x;
        
        // Area = 0.5 * ||cross||, so area^2 = 0.25 * ||cross||^2
        double crossLenSq = cross_x * cross_x + cross_y * cross_y + cross_z * cross_z;
        return crossLenSq * 0.25;
    };

    auto getLinearRootNaive = [](const Point3D& v1, const Point3D& v2, double val1, double val2) -> Point3D {
        // Avoid division by zero - if values are same, return midpoint
        if (std::abs(val1 - val2) < 1e-12) {
            return Point3D((v1.x + v2.x) * 0.5, (v1.y + v2.y) * 0.5, (v1.z + v2.z) * 0.5);
        }
        
        // Linear interpolation: t such that val1 + t*(val2 - val1) = 0
        double t = -val1 / (val2 - val1);
        
        
        return Point3D(
            v1.x + t * (v2.x - v1.x),
            v1.y + t * (v2.y - v1.y),
            v1.z + t * (v2.z - v1.z)
        );
    };
    
    // Threshold for near-degenerate triangles (squared area)
    const double minTriangleAreaSquared = 1e-16; // Very small threshold
    
    // Main marching tetrahedra function
    // Returns: {zero_crossing_points, triangles}
    auto marchTet = [&](const std::vector<Point3D>& vts, 
                        const std::vector<TetElement>& tets,
                        const std::vector<double>& vals) -> std::pair<std::vector<Point3D>, std::vector<TriangleMesh>> {
        
        std::vector<Point3D> zeroVerts;
        std::vector<TriangleMesh> triangles;
        
        // Hash map for zero-crossing points on edges (to avoid duplicates)
        // Key: sorted edge vertex indices, Value: index into zeroVerts
        std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeHash;
        edgeHash.reserve(tets.size() * 3); // Rough estimate
        
        int zeroVertCounter = 0;
        
        // Process each tet
        for (const auto& tet : tets) {
            // Get vertex indices (0-based)
            int vIndices[4] = {tet.v1, tet.v2, tet.v3, tet.v4};
            
            // Get signs at vertices
            int signs[4];
            for (int i = 0; i < 4; ++i) {
                signs[i] = getSign(vals[vIndices[i]]);
            }
            
            // Compute pattern index (binary: FromDigits reads left-to-right as MSB to LSB)
            // So pattern = signs[0]*8 + signs[1]*4 + signs[2]*2 + signs[3]*1
            int patternIdx = signs[0] * 8 + signs[1] * 4 + signs[2] * 2 + signs[3] * 1;
            
            // Get lookup table entry
            const auto& entry = tetTableEinds[patternIdx];
            
            // Skip if no zero-crossing (empty entry)
            if (entry.edgeIndices.empty()) {
                continue;
            }
            
            // Check tet orientation
            int ori = getTetOrientation(
                vts[vIndices[0]], vts[vIndices[1]], 
                vts[vIndices[2]], vts[vIndices[3]]
            );
            
            // Create zero-crossing points for each edge in the entry
            // Store local edge index -> zero vertex index mapping
            std::vector<int> localEdgeToZeroVert(6, -1);
            
            for (int edgeIdx1Based : entry.edgeIndices) {
                // Convert from 1-based to 0-based edge index
                int edgeIdx = edgeIdx1Based - 1;
                
                // Get vertex pair for this edge
                const auto& edgePair = tetEdges[edgeIdx];
                int vIdx1 = vIndices[edgePair.first];
                int vIdx2 = vIndices[edgePair.second];
                
                // Verify that signs differ (should always be true for edges in the lookup table)
                int sign1 = getSign(vals[vIdx1]);
                int sign2 = getSign(vals[vIdx2]);
                if (sign1 == sign2) {
                    // This shouldn't happen if lookup table is correct, but skip to be safe
                    continue;
                }
                
                // Create sorted edge key for hashing
                EdgeKey edgeKey(vIdx1, vIdx2);
                
                // Check if we already have a zero point for this edge
                auto it = edgeHash.find(edgeKey);
                if (it != edgeHash.end()) {
                    localEdgeToZeroVert[edgeIdx] = it->second;
                } else {
                    // Compute zero-crossing point
                    Point3D zeroPt = getLinearRootNaive(
                        vts[vIdx1], vts[vIdx2],
                        vals[vIdx1], vals[vIdx2]
                    );
                    
                    // Add to zero vertices
                    int zeroVertIdx = zeroVertCounter++;
                    zeroVerts.push_back(zeroPt);
                    edgeHash[edgeKey] = zeroVertIdx;
                    localEdgeToZeroVert[edgeIdx] = zeroVertIdx;
                }
            }
            
            // Create triangles from cycles
            for (const auto& cycle : entry.triangles) {
                // Get zero vertex indices for this triangle
                std::vector<int> triVerts;
                triVerts.reserve(3);
                
                for (int edgeIdx1Based : cycle) {
                    int edgeIdx = edgeIdx1Based - 1;
                    int zeroVertIdx = localEdgeToZeroVert[edgeIdx];
                    if (zeroVertIdx >= 0) {
                        triVerts.push_back(zeroVertIdx);
                    }
                }
                
                // Add triangle (orient based on tet orientation)
                // In Mathematica: If[ori == 1, edvtinds[[#]], Reverse[edvtinds[[#]]]]
                if (triVerts.size() == 3) {
                    if (ori == 1) {
                        // Positive orientation: use cycle order as-is
                        triangles.emplace_back(triVerts[0], triVerts[1], triVerts[2]);
                    } else {
                        // Negative or zero orientation: reverse the cycle order
                        triangles.emplace_back(triVerts[2], triVerts[1], triVerts[0]);
                    }
                }
            }
        }
        
        return {zeroVerts, triangles};
    };
    
    // Marching tet with midpoint interpolation for new refined vertices
    auto marchTetMid = [&](const std::vector<Point3D>& vts, 
                          const std::vector<TetElement>& tets,
                          const std::vector<double>& vals,
                          int originalVertexCount,
                          const std::vector<bool>& changedFromPositiveToNegative) -> std::pair<std::vector<Point3D>, std::vector<TriangleMesh>> {
        
        std::vector<Point3D> zeroVerts;
        std::vector<TriangleMesh> triangles;
        
        // Hash map for zero-crossing points on edges (to avoid duplicates)
        // Key: sorted edge vertex indices, Value: index into zeroVerts
        std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeHash;
        edgeHash.reserve(tets.size() * 3); // Rough estimate

        int zeroVertCounter = 0;
        
        // Process each tet
        for (const auto& tet : tets) {
            // Get vertex indices (0-based)
            int vIndices[4] = {tet.v1, tet.v2, tet.v3, tet.v4};
            
            // Get signs at vertices
            int signs[4];
            for (int i = 0; i < 4; ++i) {
                signs[i] = getSign(vals[vIndices[i]]);
            }
            
            // Compute pattern index (binary: FromDigits reads left-to-right as MSB to LSB)
            // So pattern = signs[0]*8 + signs[1]*4 + signs[2]*2 + signs[3]*1
            int patternIdx = signs[0] * 8 + signs[1] * 4 + signs[2] * 2 + signs[3] * 1;
            
            // Get lookup table entry
            const auto& entry = tetTableEinds[patternIdx];
            
            // Skip if no zero-crossing (empty entry)
            if (entry.edgeIndices.empty()) {
                continue;
            }
            
            // Check tet orientation
            int ori = getTetOrientation(
                vts[vIndices[0]], vts[vIndices[1]], 
                vts[vIndices[2]], vts[vIndices[3]]
            );
            
            // Create zero-crossing points for each edge in the entry
            // Store local edge index -> zero vertex index mapping
            std::vector<int> localEdgeToZeroVert(6, -1);
            
            for (int edgeIdx1Based : entry.edgeIndices) {
                // Convert from 1-based to 0-based edge index
                int edgeIdx = edgeIdx1Based - 1;
                
                // Get vertex pair for this edge
                const auto& edgePair = tetEdges[edgeIdx];
                int vIdx1 = vIndices[edgePair.first];
                int vIdx2 = vIndices[edgePair.second];
                
                // Verify that signs differ (should always be true for edges in the lookup table)
                int sign1 = getSign(vals[vIdx1]);
                int sign2 = getSign(vals[vIdx2]);
                if (sign1 == sign2) {
                    // This shouldn't happen if lookup table is correct, but skip to be safe
                    continue;
                }
                
                // Create sorted edge key for hashing
                EdgeKey edgeKey(vIdx1, vIdx2);
                
                // Check if we already have a zero point for this edge
                auto it = edgeHash.find(edgeKey);
                if (it != edgeHash.end()) {
                    localEdgeToZeroVert[edgeIdx] = it->second;
                } else {
                    // Compute zero-crossing point
                    Point3D zeroPt;
                    // If either vertex is a new refined vertex (index >= originalVertexCount),
                    // use midpoint instead of linear interpolation
                    if (((vIdx1 >= originalVertexCount) || (vIdx2 >= originalVertexCount))||(vIdx1 <= changedFromPositiveToNegative.size()&&changedFromPositiveToNegative[vIdx1]) || (vIdx2 <= changedFromPositiveToNegative.size()&&changedFromPositiveToNegative[vIdx2])) {
                        // Use midpoint
                        /*
                        zeroPt = Point3D(
                            (vts[vIdx1].x + vts[vIdx2].x) * 0.5,
                            (vts[vIdx1].y + vts[vIdx2].y) * 0.5,
                            (vts[vIdx1].z + vts[vIdx2].z) * 0.5
                        );
                        */
                        zeroPt = getLinearRoot(
                            vts[vIdx1], vts[vIdx2],
                            vals[vIdx1], vals[vIdx2]
                        );

                    } else {
                        // Use linear interpolation
                        
                        zeroPt = getLinearRoot(
                            vts[vIdx1], vts[vIdx2],
                            vals[vIdx1], vals[vIdx2]
                        );
                    }
                    
                    // Add to zero vertices
                    int zeroVertIdx = zeroVertCounter++;
                    zeroVerts.push_back(zeroPt);
                    edgeHash[edgeKey] = zeroVertIdx;
                    localEdgeToZeroVert[edgeIdx] = zeroVertIdx;
                }
            }
            
            // Create triangles from cycles
            for (const auto& cycle : entry.triangles) {
                // Get zero vertex indices for this triangle
                std::vector<int> triVerts;
                triVerts.reserve(3);
                
                for (int edgeIdx1Based : cycle) {
                    int edgeIdx = edgeIdx1Based - 1;
                    int zeroVertIdx = localEdgeToZeroVert[edgeIdx];
                    if (zeroVertIdx >= 0) {
                        triVerts.push_back(zeroVertIdx);
                    }
                }
                
                // Add triangle (orient based on tet orientation)
                // In Mathematica: If[ori == 1, edvtinds[[#]], Reverse[edvtinds[[#]]]]
                if (triVerts.size() == 3) {
                    if (ori == 1) {
                        // Positive orientation: use cycle order as-is
                        triangles.emplace_back(triVerts[0], triVerts[1], triVerts[2]);
                    } else {
                        // Negative or zero orientation: reverse the cycle order
                        triangles.emplace_back(triVerts[2], triVerts[1], triVerts[0]);
                    }
                }
            }
        }
        
        return {zeroVerts, triangles};
    };
    
    /*
    auto marchTetStartTime = std::chrono::steady_clock::now();
    // Run marching tetrahedra on original mesh (before topoMinCut)
    std::cout << "Running marching tetrahedra (before topoMinCut)..." << std::endl;
    auto [marchVerts, marchTris] = marchTet(vertices, tetrahedra, vertexAlphas);
    std::cout << "Marching tetrahedra produced " << marchVerts.size() 
              << " zero-crossing vertices and " << marchTris.size() << " triangles" << std::endl;
    
    // Report manifoldness of the marching tet result
    reportMeshManifoldness(marchVerts, marchTris, "Marching Tetrahedra (before topoMinCut)");
    
    // Write marching tet result to file
    {
        std::filesystem::create_directories("output");
        std::ofstream ofs("output/marching_tet_before_topoMinCut.obj");
        ofs << std::fixed << std::setprecision(15);
        
        // Write vertices
        for (const auto& v : marchVerts) {
            ofs << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        
        // Write triangles (OBJ uses 1-based indexing)
        for (const auto& t : marchTris) {
            ofs << "f " << (t.v1 + 1) << " " << (t.v2 + 1) << " " << (t.v3 + 1) << "\n";
        }
        
        ofs.close();
    }
    std::cout << "Wrote marching tetrahedra result to output/marching_tet_before_topoMinCut.obj" << std::endl;
    auto marchTetEndTime = std::chrono::steady_clock::now();
    auto marchTetDuration = std::chrono::duration_cast<std::chrono::milliseconds>(marchTetEndTime - marchTetStartTime).count();
    std::cout << "Marching tetrahedra took " << marchTetDuration << " ms" << std::endl;

    */
    // ============================================================================
    // END MARCHING TETRAHEDRA IMPLEMENTATION
    // ============================================================================

    // Build edge and triangle maps from tetrahedra (single pass — previously built maps twice).
    std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeMap;
    std::unordered_map<TriangleKey, int, TriangleKeyHash> triangleMap;

    const int totalVertices = static_cast<int>(vertices.size());

    {
        auto mapStartTime = std::chrono::steady_clock::now();
        const size_t nt = tetrahedra.size();
        // Upper-bound style reservation: cuts rehash churn on large meshes.
        edgeMap.reserve(std::max<size_t>(1024, nt * 3));
        triangleMap.reserve(std::max<size_t>(1024, nt * 2));
        int edgeId = totalVertices;
        for (const auto& tetElem : tetrahedra) {
            const EdgeKey edges[6] = {
                EdgeKey(tetElem.v1, tetElem.v2),
                EdgeKey(tetElem.v1, tetElem.v3),
                EdgeKey(tetElem.v1, tetElem.v4),
                EdgeKey(tetElem.v2, tetElem.v3),
                EdgeKey(tetElem.v2, tetElem.v4),
                EdgeKey(tetElem.v3, tetElem.v4)
            };
            for (int i = 0; i < 6; ++i) {
                if (edgeMap.find(edges[i]) == edgeMap.end()) {
                    edgeMap[edges[i]] = edgeId++;
                }
            }
        }
        const int totalEdgesBuilt = static_cast<int>(edgeMap.size());
        int triangleId = totalVertices + totalEdgesBuilt;
        for (const auto& tetElem : tetrahedra) {
            const TriangleKey triangles[4] = {
                TriangleKey(tetElem.v1, tetElem.v2, tetElem.v3),
                TriangleKey(tetElem.v1, tetElem.v2, tetElem.v4),
                TriangleKey(tetElem.v1, tetElem.v3, tetElem.v4),
                TriangleKey(tetElem.v2, tetElem.v3, tetElem.v4)
            };
            for (int i = 0; i < 4; ++i) {
                if (triangleMap.find(triangles[i]) == triangleMap.end()) {
                    triangleMap[triangles[i]] = triangleId++;
                }
            }
        }
        auto mapEndTime = std::chrono::steady_clock::now();
        auto mapDuration = std::chrono::duration_cast<std::chrono::milliseconds>(mapEndTime - mapStartTime).count();
        std::cout << "  Building edge and triangle maps: " << mapDuration << " ms" << std::endl;
    }

    const int totalEdges = static_cast<int>(edgeMap.size());
    const int totalTriangles = static_cast<int>(triangleMap.size());
    const int totalTets = static_cast<int>(tetrahedra.size());

    const int totalCells = totalVertices + totalEdges + totalTriangles + totalTets;

    std::vector <bool> cellInRange(totalCells, false);

    std::vector <double> cellAlpha(totalCells, 0.0);

    std::vector <int> cellDimension(totalCells, 0);
    //fill cellDimension with 0, 1, 2, 3 for vertices, edges, triangles, and tets respectively
    for (int i = 0; i < totalVertices; ++i) {
        cellDimension[i] = 0;
    }
    for (int i = totalVertices; i < totalVertices + totalEdges; ++i) {
        cellDimension[i] = 1;
    }
    for (int i = totalVertices + totalEdges; i < totalVertices + totalEdges + totalTriangles; ++i) {
        cellDimension[i] = 2;
    }
    for (int i = totalVertices + totalEdges + totalTriangles; i < totalCells; ++i) {
        cellDimension[i] = 3;
    }


    std::cout << "Assigning alpha values to vertices..." << std::endl;

    
    auto start_time = std::chrono::high_resolution_clock::now();

    
    // Assign alpha values directly from the input file
    for (int i = 0; i < totalVertices; ++i) {
        cellAlpha[i] = vertexAlphas[i];
    }
    
    std::cout << "Finished assigning alpha values to vertices" << std::endl;

    auto sectionStartTime = std::chrono::steady_clock::now();

    // Compute alpha values for edges (max of endpoint vertices)
    {
        auto startTime = std::chrono::steady_clock::now();
        for (const auto& pair : edgeMap) {
            int edge_id = pair.second;
            int v1 = pair.first.v1;
            int v2 = pair.first.v2;
            
            double edge_val = std::max(cellAlpha[v1], cellAlpha[v2]);
            cellAlpha[edge_id] = edge_val;
        }
        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Computing alpha values for edges: " << duration << " ms" << std::endl;
    }

    // Compute alpha values for triangles (max of 3 vertices)
    {
        auto startTime = std::chrono::steady_clock::now();
        for (const auto& pair : triangleMap) {
            int triangle_id = pair.second;
            int v1 = pair.first.v1;
            int v2 = pair.first.v2;
            int v3 = pair.first.v3;
            
            double triangle_val = std::max({cellAlpha[v1], cellAlpha[v2], cellAlpha[v3]});
            cellAlpha[triangle_id] = triangle_val;
        }
        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Computing alpha values for triangles: " << duration << " ms" << std::endl;
    }

    // Compute alpha values for tets (max of 4 vertices)
    {
        auto startTime = std::chrono::steady_clock::now();
        int tet_id = totalVertices + totalEdges + totalTriangles;
        for (const auto& tetElem : tetrahedra) {
            double tet_val = std::max({cellAlpha[tetElem.v1], cellAlpha[tetElem.v2], cellAlpha[tetElem.v3], cellAlpha[tetElem.v4]});
            cellAlpha[tet_id] = tet_val;
            
            tet_id++;
        }
        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Computing alpha values for tets: " << duration << " ms" << std::endl;
    }

    std::cout << "Finished setting up alpha values for all cells" << std::endl;

    
    // Now perform interface detection in the correct order
    {
        auto startTime = std::chrono::steady_clock::now();
        std::cout << "Performing inRange and interface detection..." << std::endl;

        for (int i = 0; i < totalCells; ++i) {

            cellInRange[i] = true;
                
        }
        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Interface detection: " << duration << " ms" << std::endl;
    }
    
    // ============================================================================
    // Pre-MinCut SoA core (flat arrays)
    // ============================================================================
    std::unique_ptr<std::array<double, 3>[]> verts_xyz;
    std::unique_ptr<double[]> verts_alpha;
    std::unique_ptr<std::array<int, 2>[]> edges_v2;
    std::unique_ptr<double[]> edges_alpha;
    std::unique_ptr<std::array<int, 3>[]> tris_e3;
    std::unique_ptr<double[]> tris_alpha;
    std::unique_ptr<std::array<int, 4>[]> tets_t4;
    std::unique_ptr<double[]> tets_alpha;
    size_t nVertsSoA = 0, nEdgesSoA = 0, nTrisSoA = 0, nTetsSoA = 0;

    // cell-id -> SoA index maps
    std::unique_ptr<int[]> cellIdToVertSoA(new int[totalCells]);
    std::unique_ptr<int[]> cellIdToEdgeSoA(new int[totalCells]);
    std::unique_ptr<int[]> cellIdToTriSoA(new int[totalCells]);
    std::unique_ptr<int[]> cellIdToTetSoA(new int[totalCells]);
    for (int i = 0; i < totalCells; ++i) {
        cellIdToVertSoA[i] = -1;
        cellIdToEdgeSoA[i] = -1;
        cellIdToTriSoA[i] = -1;
        cellIdToTetSoA[i] = -1;
    }

    // Build SoA arrays directly from mesh + maps (pre-MinCut canonical representation)
    {
        auto startTime = std::chrono::steady_clock::now();

        size_t nVertsIn = 0;
        for (int vid = 0; vid < totalVertices; ++vid) {
            if (cellInRange[vid]) ++nVertsIn;
        }
        size_t nEdgesIn = 0;
        for (const auto& pair : edgeMap) {
            if (cellInRange[pair.second]) ++nEdgesIn;
        }
        size_t nTrisIn = 0;
        for (const auto& pair : triangleMap) {
            if (cellInRange[pair.second]) ++nTrisIn;
        }
        size_t nTetsIn = 0;
        int soaTetCellId = totalVertices + totalEdges + totalTriangles;
        for (size_t i = 0; i < tetrahedra.size(); ++i, ++soaTetCellId) {
            if (cellInRange[soaTetCellId]) ++nTetsIn;
        }

        verts_xyz.reset(new std::array<double, 3>[nVertsIn]);
        verts_alpha.reset(new double[nVertsIn]);
        edges_v2.reset(new std::array<int, 2>[nEdgesIn]);
        edges_alpha.reset(new double[nEdgesIn]);
        tris_e3.reset(new std::array<int, 3>[nTrisIn]);
        tris_alpha.reset(new double[nTrisIn]);
        tets_t4.reset(new std::array<int, 4>[nTetsIn]);
        tets_alpha.reset(new double[nTetsIn]);
        nVertsSoA = nVertsIn;
        nEdgesSoA = nEdgesIn;
        nTrisSoA = nTrisIn;
        nTetsSoA = nTetsIn;

        size_t vWrite = 0;
        for (int vid = 0; vid < totalVertices; ++vid) {
            if (!cellInRange[vid]) continue;
            verts_xyz[vWrite] = {vertices[vid].x, vertices[vid].y, vertices[vid].z};
            verts_alpha[vWrite] = cellAlpha[vid];
            cellIdToVertSoA[vid] = static_cast<int>(vWrite);
            ++vWrite;
        }

        size_t eWrite = 0;
        for (const auto& pair : edgeMap) {
            const int edgeId = pair.second;
            if (!cellInRange[edgeId]) continue;
            const int v1 = pair.first.v1;
            const int v2 = pair.first.v2;
            const int sv1 = (v1 >= 0 && v1 < totalCells) ? cellIdToVertSoA[v1] : -1;
            const int sv2 = (v2 >= 0 && v2 < totalCells) ? cellIdToVertSoA[v2] : -1;
            if (sv1 < 0 || sv2 < 0) continue;
            edges_v2[eWrite] = {sv1, sv2};
            edges_alpha[eWrite] = cellAlpha[edgeId];
            cellIdToEdgeSoA[edgeId] = static_cast<int>(eWrite);
            ++eWrite;
        }
        nEdgesSoA = eWrite;

        size_t tWrite = 0;
        for (const auto& pair : triangleMap) {
            const int triId = pair.second;
            if (!cellInRange[triId]) continue;
            const int v1 = pair.first.v1;
            const int v2 = pair.first.v2;
            const int v3 = pair.first.v3;
            auto e1 = edgeMap.find(EdgeKey(v1, v2));
            auto e2 = edgeMap.find(EdgeKey(v1, v3));
            auto e3 = edgeMap.find(EdgeKey(v2, v3));
            if (e1 == edgeMap.end() || e2 == edgeMap.end() || e3 == edgeMap.end()) continue;
            const int se1 = cellIdToEdgeSoA[e1->second];
            const int se2 = cellIdToEdgeSoA[e2->second];
            const int se3 = cellIdToEdgeSoA[e3->second];
            if (se1 < 0 || se2 < 0 || se3 < 0) continue;
            tris_e3[tWrite] = {se1, se2, se3};
            tris_alpha[tWrite] = cellAlpha[triId];
            cellIdToTriSoA[triId] = static_cast<int>(tWrite);
            ++tWrite;
        }
        nTrisSoA = tWrite;

        size_t tetWrite = 0;
        int tetCellId = totalVertices + totalEdges + totalTriangles;
        for (const auto& tetElem : tetrahedra) {
            if (!cellInRange[tetCellId]) {
                ++tetCellId;
                continue;
            }
            auto tt1 = triangleMap.find(TriangleKey(tetElem.v1, tetElem.v2, tetElem.v3));
            auto tt2 = triangleMap.find(TriangleKey(tetElem.v1, tetElem.v2, tetElem.v4));
            auto tt3 = triangleMap.find(TriangleKey(tetElem.v1, tetElem.v3, tetElem.v4));
            auto tt4 = triangleMap.find(TriangleKey(tetElem.v2, tetElem.v3, tetElem.v4));
            if (tt1 == triangleMap.end() || tt2 == triangleMap.end() || tt3 == triangleMap.end() || tt4 == triangleMap.end()) {
                ++tetCellId;
                continue;
            }
            const int st1 = cellIdToTriSoA[tt1->second];
            const int st2 = cellIdToTriSoA[tt2->second];
            const int st3 = cellIdToTriSoA[tt3->second];
            const int st4 = cellIdToTriSoA[tt4->second];
            if (st1 < 0 || st2 < 0 || st3 < 0 || st4 < 0) {
                ++tetCellId;
                continue;
            }
            tets_t4[tetWrite] = {st1, st2, st3, st4};
            tets_alpha[tetWrite] = cellAlpha[tetCellId];
            cellIdToTetSoA[tetCellId] = static_cast<int>(tetWrite);
            ++tetWrite;
            ++tetCellId;
        }
        nTetsSoA = tetWrite;

        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Built pre-MinCut SoA core: " << duration << " ms"
                  << " (" << nVertsSoA << " v, " << nEdgesSoA << " e, "
                  << nTrisSoA << " t, " << nTetsSoA << " T)" << std::endl;
    }

    // Declare variables that will be used across multiple sections
    std::vector<Vertex> verts;
    std::vector<Edge*> edges;
    std::vector<Triangle*> triangles;
    std::vector<Tet> tets;
    std::vector<Simplex*> allSimplices;
    std::vector<int> sortedIndices;
    int totalSimplices = 0;


    // Pre-subdivision object graph build removed: pre-MinCut now uses SoA-backed cc/vals/currentVerts.
    // Keep containers declared for downstream stages that still use Simplex objects.
    
    auto sectionEndTime = std::chrono::steady_clock::now();
    auto sectionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(sectionEndTime - sectionStartTime).count();
    std::cout << "Total time for cell complex creation section: " << sectionDuration << " ms" << std::endl;

    // Cell-complex state used by refinement / TopoMinCut path.
    std::vector<std::vector<std::vector<int>>> cc(4);
    std::vector<std::vector<double>> vals(4);
    std::vector<Point3D> currentVerts;
    std::vector<std::vector<std::vector<std::vector<int>>>> finalSubs(3);

    currentVerts.resize(nVertsSoA);
    vals[0].resize(nVertsSoA);
    for (size_t i = 0; i < nVertsSoA; ++i) {
        currentVerts[i] = Point3D(verts_xyz[i][0], verts_xyz[i][1], verts_xyz[i][2]);
        vals[0][i] = verts_alpha[i];
    }

    cc[1].resize(nEdgesSoA);
    vals[1].resize(nEdgesSoA);
    for (size_t i = 0; i < nEdgesSoA; ++i) {
        cc[1][i] = {edges_v2[i][0], edges_v2[i][1]};
        vals[1][i] = edges_alpha[i];
    }

    cc[2].resize(nTrisSoA);
    vals[2].resize(nTrisSoA);
    for (size_t i = 0; i < nTrisSoA; ++i) {
        cc[2][i] = {tris_e3[i][0], tris_e3[i][1], tris_e3[i][2]};
        vals[2][i] = tris_alpha[i];
    }

    cc[3].resize(nTetsSoA);
    vals[3].resize(nTetsSoA);
    for (size_t i = 0; i < nTetsSoA; ++i) {
        cc[3][i] = {tets_t4[i][0], tets_t4[i][1], tets_t4[i][2], tets_t4[i][3]};
        vals[3][i] = tets_alpha[i];
    }

    int n = static_cast<int>(currentVerts.size());
    int k = 4;
    size_t tot = 0;
    for (int i = 1; i < k; ++i) {
        tot += cc[static_cast<size_t>(i)].size();
    }
    std::vector<Point3D> ncoords(tot);
    std::vector<double> nvals(tot, 0.0);
    int ct = 0;

    std::cout << "Initial cell complex (pre-subdivision): " << currentVerts.size() << " vertices, "
              << cc[1].size() << " edges, " << cc[2].size() << " triangles, "
              << cc[3].size() << " tets" << std::endl;

    // ============================================================================
    // INPUT TET COMPLEX SUBDIVISION ALGORITHM (+ optional subdivided simplex rebuild)
    // Toggle: comment out the next block to skip subdivision and mesh rebuild from finalSubs.
    // ============================================================================
    std::cout << "Starting input tet complex subdivision..." << std::endl;
    auto subdivisionStartTime = std::chrono::steady_clock::now();
    
    // Step 2: Identify element types: -1 = inside, 0 = interface, 1 = outside
    // Array-based classification from cc (tri->tet incidence, edge->tri incidence, vertex->edge incidence).
    std::vector<int> faceTypes(cc[2].size(), 2);   // triangles
    std::vector<int> edgeTypes(cc[1].size(), 2);   // edges
    std::vector<int> vertexTypes(currentVerts.size(), 2);
    auto getTetLabelByIndex = [&](size_t tetIdx) -> int {
        return (tetIdx < tetLabels.size()) ? tetLabels[tetIdx] : 0;
    };
    auto startTime = std::chrono::steady_clock::now();

    const size_t triCountAll = cc[2].size();
    const size_t edgeCountAll = cc[1].size();
    const size_t vertCountAll = currentVerts.size();
    std::vector<std::array<int, 2>> triIncTet(triCountAll, {-1, -1});
    std::vector<unsigned char> triIncCnt(triCountAll, 0);
    for (size_t tetIdx = 0; tetIdx < cc[3].size(); ++tetIdx) {
        for (int triIdx : cc[3][tetIdx]) {
            if (triIdx < 0 || static_cast<size_t>(triIdx) >= triCountAll) continue;
            auto& cnt = triIncCnt[static_cast<size_t>(triIdx)];
            auto& inc = triIncTet[static_cast<size_t>(triIdx)];
            if (cnt == 0) {
                inc[0] = static_cast<int>(tetIdx);
                cnt = 1;
            } else if (cnt == 1) {
                inc[1] = static_cast<int>(tetIdx);
                cnt = 2;
            } else {
                cnt = 3;
            }
        }
    }

    // Cache triangle vertex triplets once; reused in tet mask pass.
    std::vector<std::array<int, 3>> triVertsCache(triCountAll, {-1, -1, -1});
    std::vector<unsigned char> triVertsValid(triCountAll, 0);
    for (size_t triIdx = 0; triIdx < triCountAll; ++triIdx) {
        int tmp[6];
        int c = 0;
        for (int eIdx : cc[2][triIdx]) {
            if (eIdx < 0 || static_cast<size_t>(eIdx) >= edgeCountAll) continue;
            const auto& ev = cc[1][static_cast<size_t>(eIdx)];
            if (ev.size() >= 2) {
                tmp[c++] = ev[0];
                tmp[c++] = ev[1];
            }
        }
        std::sort(tmp, tmp + c);
        c = static_cast<int>(std::unique(tmp, tmp + c) - tmp);
        if (c == 3) {
            triVertsCache[triIdx] = {tmp[0], tmp[1], tmp[2]};
            triVertsValid[triIdx] = 1;
        }
    }

    for (size_t triIdx = 0; triIdx < cc[2].size(); ++triIdx) {
        const unsigned char cnt = triIncCnt[triIdx];
        const auto& inc = triIncTet[triIdx];
        if (cnt == 1) {
            faceTypes[triIdx] = 1;
        } else if (cnt == 2) {
            bool hasInside = false, hasOutside = false;
            for (int tIdx : {inc[0], inc[1]}) {
                const int label = getTetLabelByIndex(static_cast<size_t>(tIdx));
                if (label != 1) hasInside = true; else hasOutside = true;
            }
            if (hasInside && hasOutside) faceTypes[triIdx] = 0;
            else if (hasInside) faceTypes[triIdx] = -1;
            else faceTypes[triIdx] = 1;
        } else {
            faceTypes[triIdx] = 1;
        }
    }
    std::vector<unsigned char> edgeOnSurface(edgeCountAll, 0);
    std::vector<unsigned char> edgeHasInside(edgeCountAll, 0);
    std::vector<unsigned char> edgeHasOutside(edgeCountAll, 0);
    for (size_t triIdx = 0; triIdx < cc[2].size(); ++triIdx) {
        const int faceType = faceTypes[triIdx];
        for (int eIdx : cc[2][triIdx]) {
            if (eIdx < 0 || static_cast<size_t>(eIdx) >= edgeCountAll) continue;
            if (faceType == 0) edgeOnSurface[static_cast<size_t>(eIdx)] = 1;
            else if (faceType == -1) edgeHasInside[static_cast<size_t>(eIdx)] = 1;
            else if (faceType == 1) edgeHasOutside[static_cast<size_t>(eIdx)] = 1;
        }
    }
    for (size_t edgeIdx = 0; edgeIdx < edgeCountAll; ++edgeIdx) {
        const bool onSurface = edgeOnSurface[edgeIdx] != 0;
        const bool hasInsideParent = edgeHasInside[edgeIdx] != 0;
        const bool hasOutsideParent = edgeHasOutside[edgeIdx] != 0;
        if (onSurface) edgeTypes[edgeIdx] = 0;
        else if (hasInsideParent && !hasOutsideParent) edgeTypes[edgeIdx] = -1;
        else if (hasOutsideParent && !hasInsideParent) edgeTypes[edgeIdx] = 1;
        else edgeTypes[edgeIdx] = 1;
    }

    std::vector<unsigned char> vertOnSurface(vertCountAll, 0);
    std::vector<unsigned char> vertHasInside(vertCountAll, 0);
    std::vector<unsigned char> vertHasOutside(vertCountAll, 0);
    for (size_t edgeIdx = 0; edgeIdx < edgeCountAll; ++edgeIdx) {
        const auto& ev = cc[1][edgeIdx];
        if (ev.size() >= 2) {
            const int et = edgeTypes[edgeIdx];
            for (int k2 = 0; k2 < 2; ++k2) {
                const int vIdx = ev[static_cast<size_t>(k2)];
                if (vIdx < 0 || static_cast<size_t>(vIdx) >= vertCountAll) continue;
                if (et == 0) vertOnSurface[static_cast<size_t>(vIdx)] = 1;
                else if (et == -1) vertHasInside[static_cast<size_t>(vIdx)] = 1;
                else if (et == 1) vertHasOutside[static_cast<size_t>(vIdx)] = 1;
            }
        }
    }
    for (size_t vertIdx = 0; vertIdx < vertCountAll; ++vertIdx) {
        const bool onSurface = vertOnSurface[vertIdx] != 0;
        const bool hasInsideParent = vertHasInside[vertIdx] != 0;
        const bool hasOutsideParent = vertHasOutside[vertIdx] != 0;
        if (onSurface) vertexTypes[vertIdx] = 0;
        else if (hasInsideParent && !hasOutsideParent) vertexTypes[vertIdx] = -1;
        else if (hasOutsideParent && !hasInsideParent) vertexTypes[vertIdx] = 1;
        else vertexTypes[vertIdx] = 1;
    }
        
        // Create masks for elements that need to be subdivided
        // Criteria: NOT surface elements AND (consist of only surface vertices OR are parent of subdivided element)
        std::vector<bool> edgeSubdivideMask(cc[1].size(), false);
        std::vector<bool> triangleSubdivideMask(cc[2].size(), false);
        std::vector<bool> tetSubdivideMask(cc[3].size(), false);
        
        
            startTime = std::chrono::steady_clock::now();
            
            // Process edges (dimension 1) - check if not surface and all vertices are surface
            for (size_t i = 0; i < cc[1].size(); ++i) {
                if (i >= edgeTypes.size() || edgeTypes[i] == 0) {
                    continue;  // Skip surface edges
                }
                
                bool allSurfaceVertices = true;
                const auto& edgeVerts = cc[1][i];
                for (int vIdx : edgeVerts) {
                    if (vIdx < 0 || static_cast<size_t>(vIdx) >= vertexTypes.size() || vertexTypes[static_cast<size_t>(vIdx)] != 0) {
                        allSurfaceVertices = false;
                        break;
                    }
                }
                
                if (allSurfaceVertices && edgeVerts.size() == 2) {
                    edgeSubdivideMask[i] = true;
                }
            }
            
            // Process triangles (dimension 2) - check if not surface and (all vertices are surface OR has subdivided child)
            for (size_t i = 0; i < cc[2].size(); ++i) {
                if (i >= faceTypes.size() || faceTypes[i] == 0) {
                    continue;  // Skip surface triangles
                }
                
                bool allSurfaceVertices = true;
                bool hasSubdividedChild = false;
                
                for (int edgeIdx : cc[2][i]) {
                    if (edgeIdx >= 0 && static_cast<size_t>(edgeIdx) < edgeSubdivideMask.size() && edgeSubdivideMask[static_cast<size_t>(edgeIdx)]) {
                        hasSubdividedChild = true;
                    }
                    if (!hasSubdividedChild && edgeIdx >= 0 && static_cast<size_t>(edgeIdx) < cc[1].size()) {
                        for (int vIdx : cc[1][static_cast<size_t>(edgeIdx)]) {
                            if (vIdx < 0 || static_cast<size_t>(vIdx) >= vertexTypes.size() || vertexTypes[static_cast<size_t>(vIdx)] != 0) {
                                allSurfaceVertices = false;
                                break;
                            }
                        }
                    }
                }
                
                if ((allSurfaceVertices || hasSubdividedChild) && cc[2][i].size() == 3) {
                    triangleSubdivideMask[i] = true;
                }
            }
            
            // Process tets (dimension 3) - check if all vertices are surface OR has subdivided child
            for (size_t i = 0; i < cc[3].size(); ++i) {
                bool allSurfaceVertices = true;
                bool hasSubdividedChild = false;
                
                for (int triIdx : cc[3][i]) {
                    if (triIdx >= 0 && static_cast<size_t>(triIdx) < triangleSubdivideMask.size() && triangleSubdivideMask[static_cast<size_t>(triIdx)]) {
                        hasSubdividedChild = true;
                    }
                    if (!hasSubdividedChild && triIdx >= 0 && static_cast<size_t>(triIdx) < cc[2].size()) {
                        if (!triVertsValid[static_cast<size_t>(triIdx)]) {
                            allSurfaceVertices = false;
                            continue;
                        }
                        const auto& triVerts = triVertsCache[static_cast<size_t>(triIdx)];
                        for (int tv = 0; tv < 3; ++tv) {
                            const int vIdx = triVerts[static_cast<size_t>(tv)];
                            if (vIdx < 0 || static_cast<size_t>(vIdx) >= vertexTypes.size() || vertexTypes[static_cast<size_t>(vIdx)] != 0) {
                                allSurfaceVertices = false;
                                break;
                            }
                        }
                    }
                }
                
                if ((allSurfaceVertices || hasSubdividedChild) && cc[3][i].size() == 4) {
                    tetSubdivideMask[i] = true;
                }
            
            
            
        }

        // Count elements to be subdivided
        size_t numEdgesToSubdivide = 0, numTrianglesToSubdivide = 0, numTetsToSubdivide = 0;
        for (bool b : edgeSubdivideMask) if (b) numEdgesToSubdivide++;
        for (bool b : triangleSubdivideMask) if (b) numTrianglesToSubdivide++;
        for (bool b : tetSubdivideMask) if (b) numTetsToSubdivide++;


        auto endTime1 = std::chrono::steady_clock::now();
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(endTime1 - startTime).count();
        std::cout << "  Created subdivision masks: " << numEdgesToSubdivide << " edges, " 
                  << numTrianglesToSubdivide << " triangles, " << numTetsToSubdivide << " tets: " 
                  << duration1 << " ms" << std::endl;
        
        // Count elements by type
        size_t numSurfaceFaces = 0, numInsideFaces = 0, numOutsideFaces = 0;
        size_t numSurfaceEdges = 0, numInsideEdges = 0, numOutsideEdges = 0;
        size_t numSurfaceVertices = 0, numInsideVertices = 0, numOutsideVertices = 0;
        
        for (int type : faceTypes) {
            if (type == 0) numSurfaceFaces++;
            else if (type == -1) numInsideFaces++;
            else if (type == 1) numOutsideFaces++;
        }
        for (int type : edgeTypes) {
            if (type == 0) numSurfaceEdges++;
            else if (type == -1) numInsideEdges++;
            else if (type == 1) numOutsideEdges++;
        }
        for (int type : vertexTypes) {
            if (type == 0) numSurfaceVertices++;
            else if (type == -1) numInsideVertices++;
            else if (type == 1) numOutsideVertices++;
        }
        
        auto endTime = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << "  Identified faces: " << numSurfaceFaces << " surface, " 
                  << numInsideFaces << " inside, " << numOutsideFaces << " outside" << std::endl;
        std::cout << "  Identified edges: " << numSurfaceEdges << " surface, " 
                  << numInsideEdges << " inside, " << numOutsideEdges << " outside" << std::endl;
        std::cout << "  Identified vertices: " << numSurfaceVertices << " surface, " 
                  << numInsideVertices << " inside, " << numOutsideVertices << " outside" << std::endl;
        std::cout << "  Classification took: " << duration << " ms" << std::endl;
    
    
    // Helper function to compute centroid
    auto computeCentroid = [](const std::vector<std::vector<double>>& verts) -> std::vector<double> {
        std::vector<double> centroid(3, 0.0);
        for (const auto& v : verts) {
            centroid[0] += v[0];
            centroid[1] += v[1];
            centroid[2] += v[2];
        }
        double n = static_cast<double>(verts.size());
        centroid[0] /= n;
        centroid[1] /= n;
        centroid[2] /= n;
        return centroid;
    };
    
    // Helper function to compute distance between two points
    auto computeDistance = [](const std::vector<double>& p1, const std::vector<double>& p2) -> double {
        double dx = p1[0] - p2[0];
        double dy = p1[1] - p2[1];
        double dz = p1[2] - p2[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    
    // Step 3: subdivide using subs pattern (cc, vals, currentVerts, finalSubs, n, k, tot, ncoords, nvals, ct
    // are initialized before this /* block */


    std::cout << "Built cell complex: " << currentVerts.size() << " vertices, "
              << cc[1].size() << " edges, " << cc[2].size() << " triangles, "
              << cc[3].size() << " tets" << std::endl;

    // Vertex merge helper for dim>=2 with reusable scratch buffers (avoids per-cell allocations).
    std::vector<int> getVindsScratch;
    getVindsScratch.reserve(128);
    auto fillVinds = [&](int dim, int cellIdx, const std::vector<std::vector<std::vector<int>>>& currentSubs, std::vector<int>& out) {
        out.clear();
        if (dim == 1) {
            if (cellIdx < 0 || cellIdx >= static_cast<int>(cc[1].size())) return;
            out = cc[1][cellIdx];
            return;
        }
        if (cellIdx < 0 || cellIdx >= static_cast<int>(cc[dim].size())) return;
        const auto& faces = cc[dim][cellIdx];
        for (int faceIdx : faces) {
            if (faceIdx < 0 || faceIdx >= static_cast<int>(currentSubs.size())) continue;
            const auto& faceSubs = currentSubs[faceIdx];
            for (const auto& sub : faceSubs) {
                out.insert(out.end(), sub.begin(), sub.end());
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
    };
    
    // Initialize subdivisions: for each original cell, store list of refined cells (as vertex indices)
    // Start with vertex subdivisions (each vertex is itself)
    std::vector<std::vector<std::vector<int>>> subs(currentVerts.size());
    for (size_t i = 0; i < currentVerts.size(); ++i) {
        subs[i] = {{static_cast<int>(i)}};
    }
    const std::vector<std::vector<std::vector<int>>>* activeSubs = &subs;
    
    // finalSubs is declared above (outside this optional block)
    // Process cells from dimension 1 to k-1 (edges -> triangles -> tets)
    for (int i = 1; i < k; ++i) {
        int dim = i;
        std::vector<std::vector<std::vector<int>>> nsubs(cc[dim].size());
        std::cout << "Processing dim=" << dim << ", cc[" << dim << "].size()=" << cc[dim].size() << std::endl;
        
        for (size_t j = 0; j < cc[dim].size(); ++j) {
            const auto& faces = cc[dim][j];  // face indices (lower dimension)
            if (j >= vals[dim].size()) {
                continue;
            }
            
            // Determine if this cell should be subdivided using the pre-computed masks
            bool shouldSubdivide = false;
            if (dim == 1) {
                // Use edgeSubdivideMask
                if (j < edgeSubdivideMask.size()) {
                    shouldSubdivide = edgeSubdivideMask[j];
                }
            } else if (dim == 2) {
                // Use triangleSubdivideMask
                if (j < triangleSubdivideMask.size()) {
                    shouldSubdivide = triangleSubdivideMask[j];
                }
            } else if (dim == 3) {
                // Use tetSubdivideMask
                if (j < tetSubdivideMask.size()) {
                    shouldSubdivide = tetSubdivideMask[j];
                }
            }
            
            if (!shouldSubdivide) {
                // No subdivision: use face subdivisions directly
                if (dim == 1) {
                    nsubs[j] = {cc[1][j]};
                } else {
                    fillVinds(dim, static_cast<int>(j), *activeSubs, getVindsScratch);
                    if (!getVindsScratch.empty()) {
                        nsubs[j] = {getVindsScratch};
                    } else {
                        nsubs[j] = {};
                    }
                }
            } else {
                // Subdivide: create new vertex at centroid and combine with face subdivisions
                fillVinds(dim, static_cast<int>(j), *activeSubs, getVindsScratch);
                const std::vector<int>& vertIndices = getVindsScratch;
                double cx = 0.0, cy = 0.0, cz = 0.0;
                int centroidCount = 0;
                for (int vIdx : vertIndices) {
                    if (vIdx >= 0 && vIdx < static_cast<int>(currentVerts.size())) {
                        cx += currentVerts[vIdx].x;
                        cy += currentVerts[vIdx].y;
                        cz += currentVerts[vIdx].z;
                        centroidCount++;
                    }
                }
                const double inv = (centroidCount > 0) ? (1.0 / static_cast<double>(centroidCount)) : 0.0;
                const double centroidX = cx * inv;
                const double centroidY = cy * inv;
                const double centroidZ = cz * inv;
                
                // Create new vertex
                int newVertIdx = n + ct;
                ncoords[ct] = Point3D(centroidX, centroidY, centroidZ);
                
                // Compute alpha value
                double minDist = std::numeric_limits<double>::max();

                // Determine sign based on element type
                double sign = 1.0;
                if (dim == 1 && j < edgeTypes.size()) {
                    sign = (edgeTypes[j] != 1) ? -1.0 : 1.0;
                } else if (dim == 2 && j < faceTypes.size()) {
                    sign = (faceTypes[j] != 1) ? -1.0 : 1.0;
                } else if (dim == 3) {
                    // For tets, use label
                    if (j < cc[3].size()) {
                        int label = getTetLabelByIndex(j);
                        sign = (label != 1) ? -1.0 : 1.0;
                    }
                }

                int count = 0;
                for (int vIdx : vertIndices) {
                    if (vIdx >= 0 && vIdx < static_cast<int>(currentVerts.size()) && vertexTypes[vIdx] ==0) {
                        double dx = centroidX - currentVerts[vIdx].x;
                        double dy = centroidY - currentVerts[vIdx].y;
                        double dz = centroidZ - currentVerts[vIdx].z;
                        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                        dist = dist + std::abs(vals[0][vIdx]);
                        if (dist < minDist) {
                            minDist = dist;
                            count+=1;
                        }
                    }
                }

                if (count == 0) {
                    std::cout << "Error: count == 0 for dim=" << dim << ", j=" << j << std::endl;
                }
                

                nvals[ct] = sign * minDist;
                ct++;
                
                // Create new cells by combining face subdivisions with new vertex (iterate subs in place — no facesubs copy)
                size_t estimatedCells = 0;
                for (int faceIdx : faces) {
                    if (faceIdx >= 0 && faceIdx < static_cast<int>(activeSubs->size())) {
                        for (const auto& sub : (*activeSubs)[faceIdx]) {
                            if (!sub.empty()) {
                                ++estimatedCells;
                            }
                        }
                    }
                }
                std::vector<std::vector<int>> newCells;
                newCells.reserve(estimatedCells);
                for (int faceIdx : faces) {
                    if (faceIdx < 0 || faceIdx >= static_cast<int>(activeSubs->size())) {
                        continue;
                    }
                    for (const auto& sub : (*activeSubs)[faceIdx]) {
                        if (sub.empty()) {
                            continue;
                        }
                        newCells.emplace_back(sub);
                        newCells.back().push_back(newVertIdx);
                    }
                }
                nsubs[j] = std::move(newCells);
            }
        }
        
        // Keep one copy in finalSubs; move into subs for next dim (was: two full copies of nsubs).
        finalSubs[dim - 1] = std::move(nsubs);
        activeSubs = &finalSubs[dim - 1];
        std::cout << "  Created " << ct << " new vertices so far" << std::endl;
    }
    
    // Update currentVerts with new vertices
    currentVerts.insert(currentVerts.end(), ncoords.begin(), ncoords.begin() + ct);
    vals[0].insert(vals[0].end(), nvals.begin(), nvals.begin() + ct);
    
    // Count refined cells
    size_t numRefinedVertices = currentVerts.size();
    size_t numRefinedEdges = 0;
    for (const auto& edgeSub : finalSubs[0]) {
        numRefinedEdges += edgeSub.size();
    }
    size_t numRefinedTriangles = 0;
    for (const auto& triSub : finalSubs[1]) {
        numRefinedTriangles += triSub.size();
    }
    size_t numRefinedTets = 0;
    for (const auto& tetSub : finalSubs[2]) {
        numRefinedTets += tetSub.size();
    }

    std::cout << "Subdivision complete:" << std::endl;
    std::cout << "  Refined vertices: " << numRefinedVertices << std::endl;
    std::cout << "  Refined edges: " << numRefinedEdges << std::endl;
    std::cout << "  Refined triangles: " << numRefinedTriangles << std::endl;
    std::cout << "  Refined tets: " << numRefinedTets << std::endl;
    
    // Store the subdivided cell complex for use in subsequent operations
    // The subdivided complex is: currentVerts (vertices), finalSubs (subdivisions), cc (structure), vals (values)
    // We'll use this for marching tet, TopoMinCut, etc.
    
    // Note: The subdivided cell complex (cc, vals, currentVerts, finalSubs) is now ready to use
    // The old Simplex structure (allSimplices, verts, edges, triangles, tets) is replaced by
    // the subdivided cell complex structure for subsequent operations
    
    auto subdivisionEndTime = std::chrono::steady_clock::now();
    auto subdivisionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(subdivisionEndTime - subdivisionStartTime).count();
    std::cout << "Input tet complex subdivision took " << subdivisionDuration << " ms" << std::endl;

    // ============================================================================
    // END OF INPUT TET COMPLEX SUBDIVISION
    // ============================================================================
    
    // The subdivided cell complex is now stored in:
    // - currentVerts: subdivided vertices (Point3D)
    // - finalSubs: subdivisions for each dimension [edges, triangles, tets]
    // - cc: cell complex structure (will be updated with subdivided cells)
    // - vals: alpha values for all cells
    
    // Rebuild the Simplex structure from the subdivided cell complex
    // We have: vertices with values, and tets (lists of 4 vertex indices)
    // This is equivalent to the initial setup - build cell complex from scratch
    
    startTime = std::chrono::steady_clock::now();
    std::cout << "Rebuilding Simplex structure from subdivided cell complex..." << std::endl;

    // Array-only rebuild: refined vertices + refined tets -> derive edges/tris + boundary.
    const int V = static_cast<int>(currentVerts.size());
    std::vector<double> vertVals = vals[0];
    if (vertVals.size() < static_cast<size_t>(V)) vertVals.resize(static_cast<size_t>(V), 0.0);

    std::vector<std::array<int, 4>> refinedTetVerts;
    size_t numRefinedTetsExpected = 0;
    for (const auto& sub : finalSubs[2]) numRefinedTetsExpected += sub.size();
    refinedTetVerts.reserve(numRefinedTetsExpected);
    for (const auto& tetSub : finalSubs[2]) {
        for (const auto& t : tetSub) {
            if (t.size() == 4) refinedTetVerts.push_back({t[0], t[1], t[2], t[3]});
        }
    }

    std::vector<EdgeKey> edgeKeysAll;
    std::vector<TriangleKey> triKeysAll;
    edgeKeysAll.reserve(refinedTetVerts.size() * 6);
    triKeysAll.reserve(refinedTetVerts.size() * 4);
    for (const auto& tv : refinedTetVerts) {
        const int v1 = tv[0], v2 = tv[1], v3 = tv[2], v4 = tv[3];
        edgeKeysAll.emplace_back(v1, v2);
        edgeKeysAll.emplace_back(v1, v3);
        edgeKeysAll.emplace_back(v1, v4);
        edgeKeysAll.emplace_back(v2, v3);
        edgeKeysAll.emplace_back(v2, v4);
        edgeKeysAll.emplace_back(v3, v4);
        triKeysAll.emplace_back(v1, v2, v3);
        triKeysAll.emplace_back(v1, v2, v4);
        triKeysAll.emplace_back(v1, v3, v4);
        triKeysAll.emplace_back(v2, v3, v4);
    }
    std::sort(edgeKeysAll.begin(), edgeKeysAll.end());
    edgeKeysAll.erase(std::unique(edgeKeysAll.begin(), edgeKeysAll.end()), edgeKeysAll.end());
    std::sort(triKeysAll.begin(), triKeysAll.end());
    triKeysAll.erase(std::unique(triKeysAll.begin(), triKeysAll.end()), triKeysAll.end());

    const int E = static_cast<int>(edgeKeysAll.size());
    const int T = static_cast<int>(triKeysAll.size());
    const int TT = static_cast<int>(refinedTetVerts.size());
    std::vector<std::array<int, 2>> edgeVertsById(static_cast<size_t>(E));
    for (int ei = 0; ei < E; ++ei) edgeVertsById[static_cast<size_t>(ei)] = {edgeKeysAll[static_cast<size_t>(ei)].v1, edgeKeysAll[static_cast<size_t>(ei)].v2};
    std::vector<std::array<int, 3>> triVertsById(static_cast<size_t>(T));
    for (int ti = 0; ti < T; ++ti) triVertsById[static_cast<size_t>(ti)] = {triKeysAll[static_cast<size_t>(ti)].v1, triKeysAll[static_cast<size_t>(ti)].v2, triKeysAll[static_cast<size_t>(ti)].v3};
    auto edgeIdOf = [&](int a, int b) -> int {
        const EdgeKey key(a, b);
        auto it = std::lower_bound(edgeKeysAll.begin(), edgeKeysAll.end(), key);
        if (it == edgeKeysAll.end() || it->v1 != key.v1 || it->v2 != key.v2) return -1;
        return static_cast<int>(it - edgeKeysAll.begin());
    };
    auto triIdOf = [&](int a, int b, int c) -> int {
        const TriangleKey key(a, b, c);
        auto it = std::lower_bound(triKeysAll.begin(), triKeysAll.end(), key);
        if (it == triKeysAll.end() || it->v1 != key.v1 || it->v2 != key.v2 || it->v3 != key.v3) return -1;
        return static_cast<int>(it - triKeysAll.begin());
    };

    std::vector<std::array<int, 3>> triEdgesById(static_cast<size_t>(T), {-1, -1, -1});
    for (int ti = 0; ti < T; ++ti) {
        const int v1 = triVertsById[static_cast<size_t>(ti)][0];
        const int v2 = triVertsById[static_cast<size_t>(ti)][1];
        const int v3 = triVertsById[static_cast<size_t>(ti)][2];
        const int e1 = edgeIdOf(v1, v2);
        const int e2 = edgeIdOf(v1, v3);
        const int e3 = edgeIdOf(v2, v3);
        if (e1 >= 0 && e2 >= 0 && e3 >= 0) {
            triEdgesById[static_cast<size_t>(ti)] = {e1, e2, e3};
        }
    }

    std::vector<std::array<int, 4>> tetTrisById(static_cast<size_t>(TT), {-1, -1, -1, -1});
    size_t skippedInvalidTetVerts = 0, skippedMissingTriangle = 0;
    for (int tti = 0; tti < TT; ++tti) {
        const int v1 = refinedTetVerts[static_cast<size_t>(tti)][0];
        const int v2 = refinedTetVerts[static_cast<size_t>(tti)][1];
        const int v3 = refinedTetVerts[static_cast<size_t>(tti)][2];
        const int v4 = refinedTetVerts[static_cast<size_t>(tti)][3];
        if (v1 < 0 || v1 >= V || v2 < 0 || v2 >= V || v3 < 0 || v3 >= V || v4 < 0 || v4 >= V) {
            ++skippedInvalidTetVerts;
            continue;
        }
        const int t1 = triIdOf(v1, v2, v3);
        const int t2 = triIdOf(v1, v2, v4);
        const int t3 = triIdOf(v1, v3, v4);
        const int t4 = triIdOf(v2, v3, v4);
        if (t1 < 0 || t2 < 0 || t3 < 0 || t4 < 0) {
            ++skippedMissingTriangle;
            continue;
        }
        tetTrisById[static_cast<size_t>(tti)] = {t1, t2, t3, t4};
    }

    std::vector<double> edgeVals(static_cast<size_t>(E), 0.0);
    for (int ei = 0; ei < E; ++ei) {
        const int v1 = edgeVertsById[static_cast<size_t>(ei)][0];
        const int v2 = edgeVertsById[static_cast<size_t>(ei)][1];
        if (v1 >= 0 && v1 < V && v2 >= 0 && v2 < V) edgeVals[static_cast<size_t>(ei)] = std::max(vertVals[static_cast<size_t>(v1)], vertVals[static_cast<size_t>(v2)]);
    }
    std::vector<double> triVals(static_cast<size_t>(T), 0.0);
    for (int ti = 0; ti < T; ++ti) {
        const int v1 = triVertsById[static_cast<size_t>(ti)][0], v2 = triVertsById[static_cast<size_t>(ti)][1], v3 = triVertsById[static_cast<size_t>(ti)][2];
        if (v1 >= 0 && v1 < V && v2 >= 0 && v2 < V && v3 >= 0 && v3 < V) {
            triVals[static_cast<size_t>(ti)] = std::max({vertVals[static_cast<size_t>(v1)], vertVals[static_cast<size_t>(v2)], vertVals[static_cast<size_t>(v3)]});
        }
    }
    std::vector<double> tetVals(static_cast<size_t>(TT), 0.0);
    for (int tti = 0; tti < TT; ++tti) {
        const auto& tv = refinedTetVerts[static_cast<size_t>(tti)];
        if (tv[0] >= 0 && tv[0] < V && tv[1] >= 0 && tv[1] < V && tv[2] >= 0 && tv[2] < V && tv[3] >= 0 && tv[3] < V) {
            tetVals[static_cast<size_t>(tti)] = std::max({vertVals[static_cast<size_t>(tv[0])], vertVals[static_cast<size_t>(tv[1])], vertVals[static_cast<size_t>(tv[2])], vertVals[static_cast<size_t>(tv[3])]});
        }
    }

    const int N = V + E + T + TT;
    std::vector<int> dimsFlat(static_cast<size_t>(N), 0);
    std::vector<double> valsFlat(static_cast<size_t>(N), 0.0);
    std::vector<std::array<int, 4>> childrenOldFixed(static_cast<size_t>(N), std::array<int, 4>{-1, -1, -1, -1});
    std::vector<unsigned char> childCountOld(static_cast<size_t>(N), 0);
    for (int i = 0; i < V; ++i) valsFlat[static_cast<size_t>(i)] = vertVals[static_cast<size_t>(i)];
    for (int ei = 0; ei < E; ++ei) {
        const int col = V + ei;
        dimsFlat[static_cast<size_t>(col)] = 1;
        valsFlat[static_cast<size_t>(col)] = edgeVals[static_cast<size_t>(ei)];
        childrenOldFixed[static_cast<size_t>(col)][0] = edgeVertsById[static_cast<size_t>(ei)][0];
        childrenOldFixed[static_cast<size_t>(col)][1] = edgeVertsById[static_cast<size_t>(ei)][1];
        childCountOld[static_cast<size_t>(col)] = 2;
    }
    for (int ti = 0; ti < T; ++ti) {
        const int col = V + E + ti;
        dimsFlat[static_cast<size_t>(col)] = 2;
        valsFlat[static_cast<size_t>(col)] = triVals[static_cast<size_t>(ti)];
        const auto& te = triEdgesById[static_cast<size_t>(ti)];
        unsigned char cnt = 0;
        if (te[0] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + te[0];
        if (te[1] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + te[1];
        if (te[2] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + te[2];
        childCountOld[static_cast<size_t>(col)] = cnt;
    }
    for (int tti = 0; tti < TT; ++tti) {
        const int col = V + E + T + tti;
        dimsFlat[static_cast<size_t>(col)] = 3;
        valsFlat[static_cast<size_t>(col)] = tetVals[static_cast<size_t>(tti)];
        const auto& tt = tetTrisById[static_cast<size_t>(tti)];
        unsigned char cnt = 0;
        if (tt[0] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + E + tt[0];
        if (tt[1] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + E + tt[1];
        if (tt[2] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + E + tt[2];
        if (tt[3] >= 0) childrenOldFixed[static_cast<size_t>(col)][cnt++] = V + E + tt[3];
        childCountOld[static_cast<size_t>(col)] = cnt;
    }

    std::vector<int> perm(static_cast<size_t>(N));
    std::iota(perm.begin(), perm.end(), 0);
#if defined(TOPOTET_HAS_STD_EXECUTION)
    std::sort(std::execution::par, perm.begin(), perm.end(), [&](int a, int b) {
#else
    std::sort(perm.begin(), perm.end(), [&](int a, int b) {
#endif
        if (valsFlat[static_cast<size_t>(a)] != valsFlat[static_cast<size_t>(b)]) return valsFlat[static_cast<size_t>(a)] < valsFlat[static_cast<size_t>(b)];
        if (dimsFlat[static_cast<size_t>(a)] != dimsFlat[static_cast<size_t>(b)]) return dimsFlat[static_cast<size_t>(a)] < dimsFlat[static_cast<size_t>(b)];
        return a < b;
    });
    std::vector<int> oldToSorted(static_cast<size_t>(N), -1);
    for (int i = 0; i < N; ++i) oldToSorted[static_cast<size_t>(perm[static_cast<size_t>(i)])] = i;
    std::vector<int> dimsSorted(static_cast<size_t>(N), 0);
    std::vector<double> valsSorted(static_cast<size_t>(N), 0.0);
    std::vector<int> childrenOffsets(static_cast<size_t>(N) + 1, 0);
    for (int newCol = 0; newCol < N; ++newCol) {
        const int oldCol = perm[static_cast<size_t>(newCol)];
        childrenOffsets[static_cast<size_t>(newCol) + 1] =
            childrenOffsets[static_cast<size_t>(newCol)] + static_cast<int>(childCountOld[static_cast<size_t>(oldCol)]);
    }
    std::vector<int> childrenIndices(static_cast<size_t>(childrenOffsets.back()), -1);
    for (int newCol = 0; newCol < N; ++newCol) {
        const int oldCol = perm[static_cast<size_t>(newCol)];
        dimsSorted[static_cast<size_t>(newCol)] = dimsFlat[static_cast<size_t>(oldCol)];
        valsSorted[static_cast<size_t>(newCol)] = valsFlat[static_cast<size_t>(oldCol)];
        const int begin = childrenOffsets[static_cast<size_t>(newCol)];
        int writeAt = begin;
        for (unsigned char kChild = 0; kChild < childCountOld[static_cast<size_t>(oldCol)]; ++kChild) {
            const int oldRow = childrenOldFixed[static_cast<size_t>(oldCol)][kChild];
            if (oldRow >= 0 && oldRow < N) {
                const int mapped = oldToSorted[static_cast<size_t>(oldRow)];
                if (mapped >= 0) childrenIndices[static_cast<size_t>(writeAt++)] = mapped;
            }
        }
        if (writeAt != childrenOffsets[static_cast<size_t>(newCol) + 1]) {
            childrenOffsets[static_cast<size_t>(newCol) + 1] = writeAt;
        }
    }
    for (size_t i = 1; i < childrenOffsets.size(); ++i) {
        if (childrenOffsets[i] < childrenOffsets[i - 1]) {
            childrenOffsets[i] = childrenOffsets[i - 1];
        }
    }
    if (!childrenIndices.empty()) {
        childrenIndices.resize(static_cast<size_t>(childrenOffsets.back()));
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto simplice_generation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "  Rebuilt arrays: " << V << " vertices, " << E << " edges, " << T << " triangles, " << TT
              << " tets (skipped invalid-vertex tets=" << skippedInvalidTetVerts
              << ", skipped missing-triangle tets=" << skippedMissingTriangle << ")" << std::endl;
    std::cout << "Simplice generation took " << simplice_generation_duration.count() << " ms" << std::endl;

    // "simplex /2" timing for JSONL export.
    simplex2TimeMs = static_cast<double>(simplice_generation_duration.count());
    

    // Current filter behavior is "keep all" in sorted order; avoid large FilterResult deep copies.
    std::vector<int> newToOld(static_cast<size_t>(N));
    std::iota(newToOld.begin(), newToOld.end(), 0);
    /*
    // Persist sorted inputs for TopoMinCut (ASCII mode expected by your pipeline)
    {
        std::ofstream f(output_file);
        f << dimsSorted.size() << std::endl;
        for (size_t col = 0; col < dimsSorted.size(); ++col) {
            f << dimsSorted[col];
            const int begin = childrenOffsets[col];
            const int end = childrenOffsets[col + 1];

            if (dimsSorted[col] != 0 && begin == end) {
                std::cout << "ERROR: Simplice " << col << " has dimension " << dimsSorted[col] << " and no children" << std::endl;
            }

            for (int at = begin; at < end; ++at) {
                f << " " << childrenIndices[static_cast<size_t>(at)];
            }
            f << std::endl;
        }
    }
    {
        std::ofstream f2(output_file2);
        f2 << valsSorted.size() << std::endl;
        for (double v : valsSorted) {
            f2 << v << std::endl;
        }
    }
    */
    topomincut::RunOutputs tet_topo_out;
    const size_t n_bm = dimsSorted.size();
    std::vector<Eigen::Triplet<int>> bm_triplets;
    bm_triplets.reserve(n_bm * 12);
    for (size_t col = 0; col < n_bm; ++col) {
        const int begin = childrenOffsets[col];
        const int end = childrenOffsets[col + 1];
        for (int at = begin; at < end; ++at) {
            bm_triplets.emplace_back(childrenIndices[static_cast<size_t>(at)], static_cast<int>(col), 1);
        }
    }
    Eigen::SparseMatrix<int> bm_sparse(static_cast<int>(n_bm), static_cast<int>(n_bm));
    bm_sparse.setFromTriplets(bm_triplets.begin(), bm_triplets.end());

    std::vector<int> bm_dims(n_bm);
    for (size_t c = 0; c < n_bm; ++c) {
        bm_dims[c] = dimsSorted[c];
    }

    topomincut::RunParams tm_params;
    tm_params.topK = topK;
    tm_params.core = core;
    tm_params.neighborhood = neighborhood;
    tm_params.cavitySkip = cavitySkip;
    tm_params.handleSkip = handleSkip;
    tm_params.componentSkip = componentSkip;
    std::cout << "Running TopoMinCut pipeline (in-process, Eigen inputs)..." << std::endl;
    const int tm_rc = topomincut::runFromEigenSparse(bm_sparse, bm_dims, valsSorted, tm_params, &tet_topo_out);
    if (tm_rc != 0) {
        std::cerr << "TopoMinCut pipeline failed with exit code " << tm_rc << std::endl;
        return 1;
    }

    // Free memory from original data structures that are no longer needed
    std::cout << "Freeing memory from original data structures..." << std::endl;
    childrenOldFixed.clear();
    childrenOldFixed.shrink_to_fit();
    childCountOld.clear();
    childCountOld.shrink_to_fit();
    childrenOffsets.clear();
    childrenOffsets.shrink_to_fit();
    childrenIndices.clear();
    childrenIndices.shrink_to_fit();
    std::vector<bool> changedFromPositiveToNegative;
    auto getKeptSimplex = [&](int) -> Simplex* { return nullptr; };
    // Note: marching cubes on original data removed since data was cleared for memory efficiency
    // If needed, re-read data or use a different approach

    
    std::cout << "Skipping marching cubes - not needed for tet mesh" << std::endl;

    /*
    // Call TopoMinCut program
    auto quoteArg = [](const std::string& s) -> std::string {
        return "\"" + s + "\"";
    };
    */
    std::vector<double> new_alpha_values = tet_topo_out.alphas_updated;

    // Update simplex values with new alpha values.
    // Matrix / alphas passed to TopoMinCut are in *kept* order (j = 0 .. k-1), same as newToOld / output_file2.
    // Do not use allSimplices.size() here: newToOld only has k entries; indexing past k-1 is UB.
    // TopoMinCut's alphas_updated.txt lines align with j, not with old column index i.
    {
        const size_t kKeep = newToOld.size();
        const size_t nAlpha = new_alpha_values.size();
        if (nAlpha < kKeep) {
            std::cerr << "Warning: alphas_updated.txt has " << nAlpha
                      << " values but filter kept " << kKeep << " simplices; only updating min count.\n";
        } else if (nAlpha > kKeep) {
            std::cerr << "Warning: alphas_updated.txt has " << nAlpha
                      << " values but filter kept " << kKeep << " simplices; ignoring extras.\n";
        }
        const size_t nUpdate = std::min(kKeep, nAlpha);
        std::vector<bool> changedFromPositiveToNegativeVals(vertVals.size(), false);
        for (size_t j = 0; j < nUpdate; ++j) {
            const int sortedIdx = newToOld[j];
            if (sortedIdx < 0 || static_cast<size_t>(sortedIdx) >= valsSorted.size()) continue;
            const double old_val = valsSorted[static_cast<size_t>(sortedIdx)];
            const double new_val = new_alpha_values[j];
            if (!sameSign(old_val, new_val)) {
                valsSorted[static_cast<size_t>(sortedIdx)] = new_val;
                const int oldId = perm[static_cast<size_t>(sortedIdx)];
                if (oldId >= 0 && oldId < V && old_val > 0.0 && new_val < 0.0) {
                    changedFromPositiveToNegativeVals[static_cast<size_t>(oldId)] = true;
                }
            }
        }
        for (size_t sidx = 0; sidx < valsSorted.size(); ++sidx) {
            const int oldId = perm[sidx];
            if (oldId < V) vertVals[static_cast<size_t>(oldId)] = valsSorted[sidx];
            else if (oldId < V + E) edgeVals[static_cast<size_t>(oldId - V)] = valsSorted[sidx];
            else if (oldId < V + E + T) triVals[static_cast<size_t>(oldId - V - E)] = valsSorted[sidx];
            else tetVals[static_cast<size_t>(oldId - V - E - T)] = valsSorted[sidx];
        }
        vals[0] = vertVals;
        vals[1] = edgeVals;
        vals[2] = triVals;
        vals[3] = tetVals;
        // Persist for marching/snap logic later.
        changedFromPositiveToNegative.assign(changedFromPositiveToNegativeVals.begin(), changedFromPositiveToNegativeVals.end());
    }
    std::cout << "Finished updating simplex values with new alpha values" << std::endl;

    // ============================================================================
    // REFINEMENT AND SNAPPING IMPLEMENTATION (after topoMinCut)
    // ============================================================================
    

    
    std::cout << "Starting refinement and snapping after topoMinCut..." << std::endl;
    auto pipelineStartTime = std::chrono::steady_clock::now();
    auto stageStartTime = pipelineStartTime;
    auto logStage = [&](const std::string& name) {
        auto stageEnd = std::chrono::steady_clock::now();
        auto stageMs = std::chrono::duration_cast<std::chrono::milliseconds>(stageEnd - stageStartTime).count();
        std::cout << "  [PROFILE] " << name << ": " << stageMs << " ms" << std::endl;
        stageStartTime = stageEnd;
    };
    
    // Array-native handoff (no object extraction).
    std::vector<Point3D> refinedVerts = currentVerts;
    std::vector<double> refinedVals = vertVals;
    std::vector<std::vector<int>> refinedEdges(static_cast<size_t>(E));
    std::vector<double> refinedEdgeVals = edgeVals;
    for (int ei = 0; ei < E; ++ei) refinedEdges[static_cast<size_t>(ei)] = {edgeVertsById[static_cast<size_t>(ei)][0], edgeVertsById[static_cast<size_t>(ei)][1]};
    std::vector<std::vector<int>> refinedTriangles(static_cast<size_t>(T));
    std::vector<double> refinedTriangleVals = triVals;
    for (int ti = 0; ti < T; ++ti) refinedTriangles[static_cast<size_t>(ti)] = {triVertsById[static_cast<size_t>(ti)][0], triVertsById[static_cast<size_t>(ti)][1], triVertsById[static_cast<size_t>(ti)][2]};
    std::vector<std::vector<int>> refinedTets(static_cast<size_t>(TT));
    std::vector<double> refinedTetVals = tetVals;
    for (int tti = 0; tti < TT; ++tti) refinedTets[static_cast<size_t>(tti)] = {refinedTetVerts[static_cast<size_t>(tti)][0], refinedTetVerts[static_cast<size_t>(tti)][1], refinedTetVerts[static_cast<size_t>(tti)][2], refinedTetVerts[static_cast<size_t>(tti)][3]};
    logStage("extract_and_reindex_cell_complex");
    
    std::cout << "Extracted cell complex: " << refinedVerts.size() << " vertices, " 
              << refinedEdges.size() << " edges, " << refinedTriangles.size() << " triangles, "
              << refinedTets.size() << " tets" << std::endl;
    
    // Build cell complex structure cc: vector of 4 vectors
    // cc[0] = vertex coordinates (Point3D)
    // cc[1] = edges, each edge represented by vertex indices
    // cc[2] = triangles, each triangle represented by edge indices
    // cc[3] = tets, each tet represented by triangle indices
    
    // Clear and reinitialize cc and vals (they were initialized in subdivision section)
    cc.clear();
    vals.clear();
    cc.resize(4);
    vals.resize(4);
    
    // cc[0] = vertex coordinates (we'll keep refinedVerts separate for coordinates)
    // vals[0] = vertex values
    vals[0] = refinedVals;
    
    // cc[1] = edges represented by vertex indices
    cc[1] = refinedEdges;
    vals[1] = refinedEdgeVals;
    std::vector<std::vector<int>> triangleFaces(refinedTriangles.size());
    for (int ti = 0; ti < T; ++ti) triangleFaces[static_cast<size_t>(ti)] = {triEdgesById[static_cast<size_t>(ti)][0], triEdgesById[static_cast<size_t>(ti)][1], triEdgesById[static_cast<size_t>(ti)][2]};
    cc[2] = triangleFaces;
    vals[2] = refinedTriangleVals;
    
    std::vector<std::vector<int>> tetFaces(refinedTets.size());
    for (int tti = 0; tti < TT; ++tti) tetFaces[static_cast<size_t>(tti)] = {tetTrisById[static_cast<size_t>(tti)][0], tetTrisById[static_cast<size_t>(tti)][1], tetTrisById[static_cast<size_t>(tti)][2], tetTrisById[static_cast<size_t>(tti)][3]};
    cc[3] = tetFaces;
    vals[3] = refinedTetVals;
    logStage("build_low_dimensional_complex");

    // Keep original cell complex for npars computation (we don't create a filtered version)
    // cc, vals, and refinedVerts remain unchanged

    // Filter out tets where all simplices have the same sign (not needed for isosurfacing)
    std::cout << "Filtering tets with uniform sign (not needed for isosurfacing)..." << std::endl;
    
    // Track which simplices are included (start with all excluded)
    std::vector<bool> vertexIncluded(refinedVerts.size(), true);
    std::vector<bool> edgeIncluded(refinedEdges.size(), true);
    std::vector<bool> triangleIncluded(refinedTriangles.size(), true);
    std::vector<bool> tetIncluded(refinedTets.size(), true);
    
    // Helper function to check if a value is positive
    auto isPositive = [](double val) { return val > 0.0; };
    
    std::vector<int> markVertex(refinedVerts.size(), -1);
    std::vector<int> markEdge(refinedEdges.size(), -1);
    std::vector<int> markTri(refinedTriangles.size(), -1);
    std::vector<int> touchedVertices;
    std::vector<int> touchedEdges;
    std::vector<int> touchedTriangles;
    touchedVertices.reserve(32);
    touchedEdges.reserve(32);
    touchedTriangles.reserve(32);

    // Go through each tet and check if all simplices have the same sign
    for (size_t tetIdx = 0; tetIdx < refinedTets.size(); ++tetIdx) {
        // Get tet value
        double tetVal = vals[3][tetIdx];
        bool tetIsPositive = isPositive(tetVal);
        
        touchedVertices.clear();
        touchedEdges.clear();
        touchedTriangles.clear();
        
        // Get vertices directly from refinedTets
        for (int vIdx : refinedTets[tetIdx]) {
            if (vIdx >= 0 && static_cast<size_t>(vIdx) < markVertex.size() && markVertex[static_cast<size_t>(vIdx)] != static_cast<int>(tetIdx)) {
                markVertex[static_cast<size_t>(vIdx)] = static_cast<int>(tetIdx);
                touchedVertices.push_back(vIdx);
            }
        }
        
        // Get triangles from cc[3]
        for (int triIdx : cc[3][tetIdx]) {
            if (triIdx < 0 || static_cast<size_t>(triIdx) >= markTri.size()) continue;
            if (markTri[static_cast<size_t>(triIdx)] != static_cast<int>(tetIdx)) {
                markTri[static_cast<size_t>(triIdx)] = static_cast<int>(tetIdx);
                touchedTriangles.push_back(triIdx);
            }
            
            // Get edges from each triangle (cc[2])
            for (int edgeIdx : cc[2][triIdx]) {
                if (edgeIdx >= 0 && static_cast<size_t>(edgeIdx) < markEdge.size() && markEdge[static_cast<size_t>(edgeIdx)] != static_cast<int>(tetIdx)) {
                    markEdge[static_cast<size_t>(edgeIdx)] = static_cast<int>(tetIdx);
                    touchedEdges.push_back(edgeIdx);
                }
            }
        }
        
        // Check if all simplices have the same sign
        bool allSameSign = true;
        
        // Check tet itself
        // (already have tetVal and tetIsPositive)
        
        // Check all vertices
        for (int vIdx : touchedVertices) {
            if (vIdx >= 0 && vIdx < static_cast<int>(vals[0].size())) {
                bool vIsPositive = isPositive(vals[0][vIdx]);
                if (vIsPositive != tetIsPositive) {
                    allSameSign = false;
                    break;
                }
            }
        }
        
        // Check all edges
        for (int edgeIdx : touchedEdges) {
            if (edgeIdx >= 0 && edgeIdx < static_cast<int>(vals[1].size())) {
                bool eIsPositive = isPositive(vals[1][edgeIdx]);
                if (eIsPositive != tetIsPositive) {
                    allSameSign = false;
                    break;
                }
            }
        }
        
        // Check all triangles
        for (int triIdx : touchedTriangles) {
            if (triIdx >= 0 && triIdx < static_cast<int>(vals[2].size())) {
                bool tIsPositive = isPositive(vals[2][triIdx]);
                if (tIsPositive != tetIsPositive) {
                    allSameSign = false;
                    break;
                }
            }
        }
        
        if (allSameSign) {
            // Include this tet and all its simplices
            tetIncluded[tetIdx] = false;
            for (int vIdx : touchedVertices) {
                vertexIncluded[vIdx] = false;
            }
            for (int edgeIdx : touchedEdges) {
                edgeIncluded[edgeIdx] = false;
            }
            for (int triIdx : touchedTriangles) {
                triangleIncluded[triIdx] = false;
            }
        
        }
        // If allSameSign is still true, this tet is excluded (already false in arrays)
    }
    logStage("uniform_sign_tet_filter");
    
    // Index maps are already set up: -1 means excluded, >= 0 means included
    // We keep the original cell complex (cc, vals, refinedVerts) and check index maps during refinement
    
    // Count how many simplices are excluded
    size_t excludedVertices = 0, excludedEdges = 0, excludedTriangles = 0, excludedTets = 0;
    for (size_t i = 0; i < vertexIncluded.size(); ++i) {
        if (!vertexIncluded[i]) excludedVertices++;
    }
    for (size_t i = 0; i < edgeIncluded.size(); ++i) {
        if (!edgeIncluded[i]) excludedEdges++;
    }
    for (size_t i = 0; i < triangleIncluded.size(); ++i) {
        if (!triangleIncluded[i]) excludedTriangles++;
    }
    for (size_t i = 0; i < tetIncluded.size(); ++i) {
        if (!tetIncluded[i]) excludedTets++;
    }
    
    std::cout << "Filtering complete: " << excludedVertices << " vertices, "
              << excludedEdges << " edges, " << excludedTriangles << " triangles, "
              << excludedTets << " tets will be skipped during refinement" << std::endl;
    
    // Keep original cell complex - we'll check index maps during refinement to skip excluded cells

    // Refinement algorithm following Mathematica pattern
    n = static_cast<int>(currentVerts.size());
    k = 4;  // dimensions: 0=vertices, 1=edges, 2=triangles, 3=tets
    
    // Allocate space for new vertex coordinates and values
    // At most one new vertex per k-cell for k>0
    tot = 0;
    for (int i = 1; i < k; ++i) {
        tot += cc[i].size();
    }

    ncoords.clear();
    ncoords.resize(tot);
    nvals.clear();
    nvals.resize(tot, 0.0);

    ct = 0;  // new vertex counter
    
    // Helper function: get vertex indices for a cell
    // For edges: directly from cc[1][j]
    // For triangles: from edge subdivisions
    // For tets: from triangle subdivisions
    // Declare as std::function so it can be reassigned

    /*
    getVinds = [&](int dim, int cellIdx, const std::vector<std::vector<std::vector<int>>>& currentSubs) -> std::vector<int> {
        if (dim == 1) {
            // Fast path for edges: directly return vertex indices
            if (cellIdx < 0 || cellIdx >= static_cast<int>(cc[1].size())) {
                std::cerr << "Error in getVinds: cellIdx " << cellIdx << " out of bounds for cc[1].size()=" << cc[1].size() << std::endl;
                return {};
            }
            return cc[1][cellIdx];  // edges are directly vertex indices
        } else {
            // For higher dimensions, get vertices from face subdivisions
            if (cellIdx < 0 || cellIdx >= static_cast<int>(cc[dim].size())) {
                std::cerr << "Error in getVinds: cellIdx " << cellIdx << " out of bounds for cc[" << dim << "].size()=" << cc[dim].size() << std::endl;
                return {};
            }
            const auto& faces = cc[dim][cellIdx];
            
            // Estimate size: for triangles (dim=2) we have 3 faces, for tets (dim=3) we have 4 faces
            // Each face typically has 1-2 subdivisions, each with 2-4 vertices
            // Reserve space to avoid rehashing
            size_t estimatedSize = (dim == 2) ? 6 : 12;  // Conservative estimate
            std::unordered_set<int> verts;
            verts.reserve(estimatedSize);
            
            for (int faceIdx : faces) {
                if (faceIdx < 0 || faceIdx >= static_cast<int>(currentSubs.size())) {
                    std::cerr << "Warning in getVinds: faceIdx " << faceIdx << " out of bounds for currentSubs.size()=" << currentSubs.size() << ", dim=" << dim << ", cellIdx=" << cellIdx << std::endl;
                    continue;
                }
                const auto& faceSubs = currentSubs[faceIdx];
                for (const auto& sub : faceSubs) {
                    // Use insert with hint for better performance when inserting multiple values
                    verts.insert(sub.begin(), sub.end());
                }
            }
            
            // Convert to vector - reserve space to avoid reallocation
            std::vector<int> result;
            result.reserve(verts.size());
            result.assign(verts.begin(), verts.end());
            return result;
        }
    };
    
    */



    auto refinementStartTime = std::chrono::steady_clock::now();

    // Initialize subdivisions: for each original cell, store list of refined cells (as vertex indices)
    // Start with vertex subdivisions (each vertex is itself)


    subs.clear();
    subs.resize(currentVerts.size());
    for (size_t i = 0; i < currentVerts.size(); ++i) {
        subs[i] = {{static_cast<int>(i)}};
    }
    const std::vector<std::vector<std::vector<int>>>* activeRefSubs = &subs;

    finalSubs.clear();
    finalSubs.resize(3);
    
    // Process cells from dimension 2 to k (edges -> triangles -> tets)
    // In Mathematica: Do[..., {i, 2, k}], so i=2 means edges (dimension 1), i=3 means triangles (dimension 2), i=4 means tets (dimension 3)
    for (int i = 1; i < k; ++i) {
        int dim = i;
        std::vector<std::vector<std::vector<int>>> nsubs(cc[dim].size());
        std::cout << "Processing dim=" << dim << ", cc[" << dim << "].size()=" << cc[dim].size() 
                  << ", subs.size()=" << activeRefSubs->size() << std::endl;
        for (size_t j = 0; j < cc[dim].size(); ++j) {
            try {
                // Early exit: Check if this cell is excluded (index map is -1)
                /*
                bool isExcluded = false;
                if (dim == 1 && j < edgeIncluded.size() && !edgeIncluded[j]) {
                    isExcluded = true;
                } else if (dim == 2 && j < triangleIncluded.size() && !triangleIncluded[j]) {
                    isExcluded = true;
                } else if (dim == 3 && j < tetIncluded.size() && !tetIncluded[j]) {
                    isExcluded = true;
                }
                
                if (isExcluded) {
                    // Cell is excluded: use original cell as its subdivision (no refinement)
                    // For edges (dim=1), directly use cc[1][j]
                    if (dim == 1) {
                        nsubs[j] = {cc[1][j]};
                    } else {
                        // For higher dimensions, need to get vertices from faces
                        std::vector<int> verts = getVinds(dim, j, subs);
                        if (!verts.empty()) {
                            nsubs[j] = {verts};
                        } else {
                            nsubs[j] = {};
                        }
                    }
                    continue;
                }
                */
                const auto& faces = cc[dim][j];  // face indices (lower dimension)
                if (j >= vals[dim].size()) {
                    std::cerr << "Error: j=" << j << " out of bounds for vals[" << dim << "].size()=" << vals[dim].size() << std::endl;
                    continue;
                }
                double val = vals[dim][j];  // cell value
                
                // Fast path: if val <= 0, cell is simple - skip expensive computations
                if (val <= 0.0) {
                    // For edges, directly use cc[1][j]
                    if (dim == 1) {
                        nsubs[j] = {cc[1][j]};
                    } else {
                        // For higher dimensions, get vertices from faces
                        fillVinds(dim, static_cast<int>(j), *activeRefSubs, getVindsScratch);
                        if (!getVindsScratch.empty()) {
                            nsubs[j] = {getVindsScratch};
                        } else {
                            nsubs[j] = {};
                        }
                    }
                    continue;
                }
                
                // Check if cell is simple: (max face value > 0 && max face subs list length == 1)
                bool isSimple = false;
                bool hasValidFaceVal = false;
                double maxFaceVal = -std::numeric_limits<double>::infinity();
                for (int faceIdx : faces) {
                    if (faceIdx < 0 || faceIdx >= static_cast<int>(vals[dim-1].size())) {
                        std::cerr << "Error: faceIdx " << faceIdx << " out of bounds for vals[" << (dim-1)
                                  << "].size()=" << vals[dim-1].size() << ", dim=" << dim << ", j=" << j << std::endl;
                        continue;
                    }
                    hasValidFaceVal = true;
                    maxFaceVal = std::max(maxFaceVal, vals[dim-1][faceIdx]);
                }
                if (hasValidFaceVal) {
                    if (maxFaceVal > 0.0) {
                        int maxSubLength = 0;
                        for (int faceIdx : faces) {
                            if (faceIdx < 0 || faceIdx >= static_cast<int>(activeRefSubs->size())) {
                                continue;
                            }
                            maxSubLength = std::max(maxSubLength, static_cast<int>((*activeRefSubs)[faceIdx].size()));
                            if (maxSubLength > 1) {
                                break;
                            }
                        }
                        isSimple = (maxSubLength == 1);
                    }
                } else {
                    std::cout << "Facevals empty for cell (this shouldnt print)" << j << std::endl;
                }
                
                // Get vertex indices only if not simple (deferred computation)
                const std::vector<int>* vertsPtr = nullptr;
                if (!isSimple) {
                    fillVinds(dim, static_cast<int>(j), *activeRefSubs, getVindsScratch);
                    if (getVindsScratch.empty()) {
                        std::cerr << "Warning: empty verts for dim=" << dim << ", j=" << j << std::endl;
                        nsubs[j] = {};
                        continue;
                    }
                    vertsPtr = &getVindsScratch;
                }
                
                if (isSimple) {
                    // Cell is simple: use original cell
                    // For edges, directly use cc[1][j]
                    if (dim == 1) {
                        nsubs[j] = {cc[1][j]};
                    } else {
                        // For higher dimensions, get vertices from faces (lazy evaluation)
                        fillVinds(dim, static_cast<int>(j), *activeRefSubs, getVindsScratch);
                        if (!getVindsScratch.empty()) {
                            nsubs[j] = {getVindsScratch};
                        } else {
                            nsubs[j] = {};
                        }
                    }
                } else {
                    // Cell is not simple: create new vertex, and form one cell connecting the vertex to each face subdivision
                    // Compute centroid (optimized: pre-compute inverse to avoid repeated division)
                    Point3D centroid(0, 0, 0);
                    const std::vector<int>& verts = *vertsPtr;
                    const size_t numVerts = verts.size();
                    const double invNumVerts = 1.0 / numVerts;
                    
                    for (int vIdx : verts) {
                        if (vIdx < static_cast<int>(refinedVerts.size())) {
                            // Original vertex
                            centroid.x += refinedVerts[vIdx].x;
                            centroid.y += refinedVerts[vIdx].y;
                            centroid.z += refinedVerts[vIdx].z;
                        } else {
                            // New refined vertex (from previous refinement steps)
                            int ncoordsIdx = vIdx - static_cast<int>(refinedVerts.size());
                            if (ncoordsIdx >= 0 && ncoordsIdx < static_cast<int>(ncoords.size())) {
                                centroid.x += ncoords[ncoordsIdx].x;
                                centroid.y += ncoords[ncoordsIdx].y;
                                centroid.z += ncoords[ncoordsIdx].z;
                            } else {
                                std::cerr << "Error: ncoordsIdx " << ncoordsIdx << " out of bounds for ncoords.size()=" 
                                          << ncoords.size() << ", vIdx=" << vIdx << ", dim=" << dim << ", j=" << j << std::endl;
                            }
                        }
                    }
                    centroid.x *= invNumVerts;
                    centroid.y *= invNumVerts;
                    centroid.z *= invNumVerts;
                    
                    if (ct >= static_cast<int>(ncoords.size())) {
                        std::cerr << "Error: ct=" << ct << " >= ncoords.size()=" << ncoords.size() 
                                  << ", dim=" << dim << ", j=" << j << std::endl;
                        throw std::runtime_error("New vertex counter out of bounds");
                    }
                    
                    ncoords[ct] = centroid;
                    nvals[ct] = val;
                    int newVertIdx = n + ct;
                    ct++;
                    
                    size_t estimatedCells = 0;
                    for (int faceIdx : faces) {
                        if (faceIdx >= 0 && faceIdx < static_cast<int>(activeRefSubs->size())) {
                            estimatedCells += (*activeRefSubs)[faceIdx].size();
                        }
                    }
                    std::vector<std::vector<int>> newCells;
                    newCells.reserve(estimatedCells);
                    for (int faceIdx : faces) {
                        if (faceIdx < 0 || faceIdx >= static_cast<int>(activeRefSubs->size())) {
                            continue;
                        }
                        for (const auto& sub : (*activeRefSubs)[faceIdx]) {
                            if (sub.empty()) {
                                std::cerr << "Warning: empty subdivision for dim=" << dim << ", j=" << j << std::endl;
                                continue;
                            }
                            newCells.emplace_back(sub);
                            newCells.back().push_back(newVertIdx);
                            if (dim == 3 && newCells.back().size() != 4) {
                                std::cerr << "Warning: tet cell has " << newCells.back().size() 
                                          << " vertices instead of 4 for dim=" << dim << ", j=" << j << std::endl;
                            }
                        }
                    }
                    nsubs[j] = std::move(newCells);
                    
                }
            } catch (const std::exception& e) {
                std::cerr << "Exception at dim=" << dim << ", j=" << j << ": " << e.what() << std::endl;
                throw;
            } catch (...) {
                std::cerr << "Unknown exception at dim=" << dim << ", j=" << j << std::endl;
                throw;
            }
        }
       
        // Store final subdivisions for this dimension (single copy: move + assign)
        // finalSubs[0] = edges, finalSubs[1] = triangles, finalSubs[2] = tets
        finalSubs[dim - 1] = std::move(nsubs);
        activeRefSubs = &finalSubs[dim - 1];
        std::cout<<"ct: "<< ct << std::endl;
    }

    // Join all subdivisions together: Apply[Join, subs] for the final dimension (tets)
    size_t refinedTetCount = 0;
    for (const auto& tetSub : *activeRefSubs) {  // final dimension subdivisions (tets)
        refinedTetCount += tetSub.size();
    }
    std::vector<std::vector<int>> refinedTetList;
    refinedTetList.resize(refinedTetCount);
    {
        size_t w = 0;
        for (const auto& tetSub : *activeRefSubs) {
            for (const auto& sub : tetSub) {
                refinedTetList[w++] = sub;
            }
        }
    }
    
    std::vector<Point3D> origVerts;
    origVerts = refinedVerts;
    // Update refinedVerts and refinedVals with new vertices
    refinedVerts.insert(refinedVerts.end(), ncoords.begin(), ncoords.begin() + ct);

    refinedVals.insert(refinedVals.end(), nvals.begin(), nvals.begin() + ct);

    size_t numRefinedVerticesAfterSnap = 0;
    size_t numRefinedEdgesAfterSnap = 0;
    size_t numRefinedTrianglesAfterSnap = 0;
    size_t numRefinedTetsAfterSnap = 0;

    // Count refined cells
    numRefinedVerticesAfterSnap = refinedVerts.size();
    numRefinedEdgesAfterSnap = 0;
    for (const auto& edgeSub : finalSubs[0]) {  // finalSubs[0] corresponds to edges (dimension 1)
        numRefinedEdgesAfterSnap += edgeSub.size();
    }
    numRefinedTrianglesAfterSnap = 0;
    for (const auto& triSub : finalSubs[1]) {  // finalSubs[1] corresponds to triangles (dimension 2)
        numRefinedTrianglesAfterSnap += triSub.size();
    }
    numRefinedTetsAfterSnap = refinedTetList.size();
    
    std::cout << "Refinement complete:" << std::endl;
    std::cout << "  Refined vertices: " << numRefinedVerticesAfterSnap << std::endl;
    std::cout << "  Refined edges: " << numRefinedEdgesAfterSnap << std::endl;
    std::cout << "  Refined triangles (faces): " << numRefinedTrianglesAfterSnap << std::endl;
    std::cout << "  Refined tets: " << numRefinedTetsAfterSnap << std::endl;
    
    // Now generate snapping mask for vertices with val <= 0
    std::vector<bool> snapMask(refinedVerts.size(), false);
    for (size_t i = 0; i < origVerts.size(); ++i) {
        if (refinedVals[i] <= 0.0) {
            snapMask[i] = true;  // Initially mark all negative vertices as potentially snappable
        }
    }

    // Count before filtering
    size_t numSnappableBefore = 0;
    for (bool b : snapMask) if (b) numSnappableBefore++;
    
    // Implement genSnapMask3D logic to filter snappable vertices
    // This filters vertices based on:
    // 1. Boundary triangles (triangles with exactly one incident tet having val <= 0)
    // 2. Connected components of boundary triangles
    // 3. Non-manifold edges (edges incident to > 2 boundary triangles)
    // 4. Isolated vertices (not part of any boundary triangle)

    // Compute npars: count how many non-positive cells of dimension i have each cell of dimension i-1 as a face
    // npars = Map[0 &, cc[[;; -2]], {2}];
    // Do[MapThread[If[#2, Map[(npars[[i - 1, #]]++) &, #1]] &, {cc[[i]], sel[[i]]}], {i, 2, k}];
    // NOTE: Use original (unfiltered) cell complex for npars computation
    int k_original = 4;  // dimensions: 0=vertices, 1=edges, 2=triangles, 3=tets
    
    // First, create sel: boolean vector for each dimension indicating if cell value <= 0
    std::vector<std::vector<bool>> sel(k_original);
    for (int dim = 0; dim < k_original; ++dim) {
        sel[dim].resize(vals[dim].size());
        for (size_t j = 0; j < vals[dim].size(); ++j) {
            sel[dim][j] = (vals[dim][j] <= 0.0);
        }
    }
    
    // Initialize npars: for dimensions 0 to k-2 (exclude last dimension)
    // npars[i] corresponds to dimension i, and npars[i][j] counts how many cells of dimension i+1 (that are <= 0) have cell j as a face
    std::vector<std::vector<int>> npars(k_original - 1);
    npars[0].resize(refinedVerts.size(), 0);
    for (int dim = 1; dim < k_original - 1; ++dim) {
        npars[dim].resize(cc[dim].size(), 0);
    }
    
    // Do loop: iterate from i=2 to k (Mathematica 1-indexed, so dimensions 1, 2, 3 in 0-indexed)
    // For each dimension i, if sel[i][j] is true, increment npars[i-1][faceIdx] for each face in cc[i][j]
    for (int i = 1; i < k_original; ++i) {  // i=1,2,3 (edges, triangles, tets)
        for (size_t j = 0; j < cc[i].size(); ++j) {
            if (sel[i][j]) {  // If cell j of dimension i has value <= 0
                // For each face index in cc[i][j], increment npars[i-1][faceIdx]
                for (int faceIdx : cc[i][j]) {
                    if (faceIdx >= 0 && faceIdx < static_cast<int>(npars[i-1].size())) {
                        npars[i-1][faceIdx]++;
                    }
                }
            }
        }
    }
    
    // Print npars for debugging
    std::cout << "npars computed:" << std::endl;
    for (int dim = 0; dim < k - 1; ++dim) {
        std::cout << "  Dimension " << dim << " (npars[" << dim << "]): ";
        size_t nonZeroCount = 0;
        for (size_t j = 0; j < npars[dim].size(); ++j) {
            if (npars[dim][j] > 0) {
                nonZeroCount++;
            }
        }
        std::cout << nonZeroCount << " non-zero entries out of " << npars[dim].size() << " total" << std::endl;
    }

    // Build vinds: vertex indices for each cell of each dimension (for dimensions 0 to k-2)
    // vinds[dim][j] contains the vertex indices for cell j of dimension dim
    std::vector<std::vector<std::vector<int>>> vinds(k - 1);
    
    // Dimension 0 (vertices): each vertex is itself
    vinds[0].resize(origVerts.size());
    for (size_t j = 0; j < origVerts.size(); ++j) {
        vinds[0][j] = {static_cast<int>(j)};
    }
    
    // Dimension 1 (edges): edges are directly vertex indices
    vinds[1].resize(cc[1].size());
    for (size_t j = 0; j < cc[1].size(); ++j) {
        vinds[1][j] = cc[1][j];
    }
    
    // Dimension 2 (triangles): get vertices from edges (fixed-size local buffer, no hash allocations)
    vinds[2].resize(cc[2].size());
    for (size_t j = 0; j < cc[2].size(); ++j) {
        int tmp[6];
        int cnt = 0;
        for (int edgeIdx : cc[2][j]) {
            if (edgeIdx >= 0 && edgeIdx < static_cast<int>(cc[1].size())) {
                const auto& ev = cc[1][edgeIdx];
                if (ev.size() >= 2) {
                    tmp[cnt++] = ev[0];
                    tmp[cnt++] = ev[1];
                }
            }
        }
        std::sort(tmp, tmp + cnt);
        cnt = static_cast<int>(std::unique(tmp, tmp + cnt) - tmp);
        auto& out = vinds[2][j];
        out.resize(static_cast<size_t>(cnt));
        for (int k2 = 0; k2 < cnt; ++k2) out[static_cast<size_t>(k2)] = tmp[k2];
    }
    
    // Mark all vertices of isolated sel-True cells non-snappable
    // MapThread[MapThread[If[#2 && #3 == 0, Map[(mask[[#]] = False) &, #1]] &, {#1, #2, #3}] &, {vinds[[;; -2]], sel[[;; -2]], npars}];
    for (int dim = 0; dim < k - 1; ++dim) {  // dimensions 0, 1, 2 (exclude last dimension)
        for (size_t j = 0; j < vinds[dim].size(); ++j) {
            // If sel[dim][j] is true (cell has value <= 0) AND npars[dim][j] == 0 (isolated)
            if (j < sel[dim].size() && sel[dim][j] && 
                j < npars[dim].size() && npars[dim][j] == 0) {
                // Mark all vertices in vinds[dim][j] as non-snappable
                for (int vIdx : vinds[dim][j]) {
                    if (vIdx >= 0 && vIdx < static_cast<int>(snapMask.size())) {
                        snapMask[vIdx] = false;
                    }
                }
            }
        }
    }
    
    std::cout << "Marked vertices of isolated sel-True cells as non-snappable" << std::endl;

    // Get non-manifold exterior edges and vertices
    // extrinds = Select[Range[Length[cc[[3]]]], (npars[[3, #]] == 1) &];
    // In Mathematica: cc[[3]] is triangles (dimension 2 in 0-indexed), npars[[3, #]] is npars[2][#] in C++
    std::vector<int> extrinds;
    extrinds.reserve(cc[2].size());
    for (size_t j = 0; j < cc[2].size(); ++j) {  // cc[2] is triangles
        if (j < npars[2].size() && npars[2][j] == 1) {
            extrinds.push_back(static_cast<int>(j));
        }
    }
    std::cout << "Found " << extrinds.size() << " exterior triangles (npars[2][j] == 1)" << std::endl;

    // Get vertex indices for exterior triangles: vinds[[3, extrinds]] in Mathematica
    // In Mathematica: vinds[[3]] is triangles (dimension 2 in 0-indexed)
    // So we select vinds[2][extrinds[i]] for each i
    std::vector<std::vector<int>> exteriorTriangles;
    exteriorTriangles.reserve(extrinds.size());
    for (int triIdx : extrinds) {
        if (triIdx >= 0 && triIdx < static_cast<int>(vinds[2].size())) {
            exteriorTriangles.push_back(vinds[2][triIdx]);
        }
    }
    
    // getSimplicialNonmaniEdges: find non-manifold edges from simplices
    // For triangles, dropping each vertex gives edges (faces of dimension d-1)
    auto getSimplicialNonmaniEdges = [](const std::vector<std::vector<int>>& simps) -> std::vector<std::vector<int>> {
        if (simps.empty()) return {};
        
        int d = static_cast<int>(simps[0].size());  // dimension of simplices (3 for triangles)
        int n = static_cast<int>(simps.size());
        
        // Hash table to map sorted face to index
        struct FaceKey {
            std::vector<int> vertices;
            bool operator==(const FaceKey& other) const {
                return vertices == other.vertices;
            }
        };
        struct FaceKeyHash {
            size_t operator()(const FaceKey& key) const {
                size_t h = 0;
                for (int v : key.vertices) {
                    h ^= std::hash<int>()(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
                }
                return h;
            }
        };
        
        std::unordered_map<FaceKey, int, FaceKeyHash> hash;
        std::vector<std::vector<int>> faces;
        faces.reserve(simps.size() * d);
        std::vector<int> degs;
        degs.reserve(simps.size() * d);
        int ct = 0;
        
        // Process each simplex
        for (const auto& simp : simps) {
            // Drop each vertex to get faces
            for (int i = 0; i < d; ++i) {
                std::vector<int> f;
                f.reserve(3);
                for (int j = 0; j < d; ++j) {
                    if (j != i) {
                        f.push_back(simp[j]);
                    }
                }
                std::sort(f.begin(), f.end());
                
                FaceKey key;
                key.vertices = f;
                
                auto it = hash.find(key);
                int ind;
                if (it == hash.end()) {
                    ct++;
                    faces.push_back(f);
                    degs.push_back(0);
                    hash[key] = ct - 1;
                    ind = ct - 1;
                } else {
                    ind = it->second;
                }
                degs[ind]++;
            }
        }
        
        // Select faces with degree != 2
        std::vector<std::vector<int>> nonManifoldEdges;
        nonManifoldEdges.reserve(ct);
        for (int i = 0; i < ct; ++i) {
            if (degs[i] != 2) {
                nonManifoldEdges.push_back(faces[i]);
            }
        }
        nonManifoldEdges.shrink_to_fit();
        
        return nonManifoldEdges;
    };
    
    // getPolygonalNonmaniVerts: find non-manifold vertices from polygons
    auto getPolygonalNonmaniVerts = [](const std::vector<std::vector<int>>& polys) -> std::vector<int> {
        if (polys.empty()) return {};
        
        // Find max vertex index
        int m = 0;
        for (const auto& poly : polys) {
            for (int v : poly) {
                m = std::max(m, v);
            }
        }
        
        // faces[i] stores edges incident to vertex i
        std::vector<std::vector<std::pair<int, int>>> faces(m + 1);
        
        // For each polygon, create edges for each vertex
        // MapThread[AppendTo[faces[[#1]], #2 \[UndirectedEdge] #3] &, {#, RotateLeft[#], RotateRight[#]}]
        // For vertex at position i, create edge between RotateLeft[i] and RotateRight[i]
        for (const auto& poly : polys) {
            int n = static_cast<int>(poly.size());
            for (int i = 0; i < n; ++i) {
                int v = poly[i];
                int vLeft = poly[(i + 1) % n];  // RotateLeft
                int vRight = poly[(i + n - 1) % n];  // RotateRight
                
                // Store edge vLeft-vRight incident to vertex v
                if (v >= 0 && v <= m) {
                    // Store as undirected edge (sorted)
                    int e1 = std::min(vLeft, vRight);
                    int e2 = std::max(vLeft, vRight);
                    faces[v].push_back({e1, e2});
                }
            }
        }
        
        // For each vertex, check if graph of edges has more than one connected component
        std::vector<int> nonManifoldVerts;
        for (int v = 0; v <= m; ++v) {
            if (faces[v].empty()) continue;
            
            // Build graph from edges
            std::unordered_map<int, std::vector<int>> graph;
            for (const auto& edge : faces[v]) {
                graph[edge.first].push_back(edge.second);
                graph[edge.second].push_back(edge.first);
            }
            
            // Find connected components using BFS
            std::unordered_set<int> visited;
            int numComponents = 0;
            
            for (const auto& [node, neighbors] : graph) {
                if (visited.find(node) != visited.end()) continue;
                
                // BFS
                std::queue<int> q;
                q.push(node);
                visited.insert(node);
                numComponents++;
                
                while (!q.empty()) {
                    int current = q.front();
                    q.pop();
                    
                    auto it = graph.find(current);
                    if (it != graph.end()) {
                        for (int neighbor : it->second) {
                            if (visited.find(neighbor) == visited.end()) {
                                visited.insert(neighbor);
                                q.push(neighbor);
                            }
                        }
                    }
                }
            }
            
            // If more than one connected component, vertex is non-manifold
            if (numComponents > 1) {
                nonManifoldVerts.push_back(v);
            }
        }
        
        return nonManifoldVerts;
    };
    
    // Compute non-manifold edges and vertices
    std::vector<std::vector<int>> nmeds = getSimplicialNonmaniEdges(exteriorTriangles);
    std::vector<int> nmvts = getPolygonalNonmaniVerts(exteriorTriangles);
    
    std::cout << "Found " << nmeds.size() << " non-manifold exterior edges" << std::endl;
    std::cout << "Found " << nmvts.size() << " non-manifold exterior vertices" << std::endl;

    // Mark non-manifold exterior edges and vertices as non-snappable
    // Map[(mask[[#]] = {False, False}) &, nmeds];  // For each edge, mark both vertices as False
    for (const auto& edge : nmeds) {
        for (int vIdx : edge) {
            if (vIdx >= 0 && vIdx < static_cast<int>(snapMask.size())) {
                snapMask[vIdx] = false;
            }
        }
    }
    
    // Map[(mask[[#]] = False) &, nmvts];  // For each vertex, mark as False
    for (int vIdx : nmvts) {
        if (vIdx >= 0 && vIdx < static_cast<int>(snapMask.size())) {
            snapMask[vIdx] = false;
        }
    }
    
    std::cout << "Marked " << nmeds.size() << " non-manifold exterior edges and " 
              << nmvts.size() << " non-manifold exterior vertices as non-snappable" << std::endl;

    
    // Count after filtering
    size_t numSnappableAfter = 0;
    for (bool b : snapMask) if (b) numSnappableAfter++;
    
    std::cout << "  Snappable vertices: " << numSnappableBefore << " -> " << numSnappableAfter << " (after filtering)" << std::endl;
    
    // Marching tet with snapping on refined mesh
    std::cout << "Running marching tetrahedra with snapping on refined mesh..." << std::endl;
    
    // Convert refined tets to TetElement format
    std::vector<TetElement> refinedTetElements;
    refinedTetElements.reserve(refinedTetList.size());
    for (const auto& tet : refinedTetList) {
        if (tet.size() == 4) {
            refinedTetElements.emplace_back(tet[0], tet[1], tet[2], tet[3]);
        }
    }
    refinedTetElements.shrink_to_fit();
    
    // Marching tet with snapping
    auto marchTetSnap = [&](const std::vector<Point3D>& vts,
                            const std::vector<TetElement>& tets,
                            const std::vector<double>& vals,
                            const std::vector<bool>& mask,
                            int originalVertexCount,
                            const std::vector<bool>& changedFromPositiveToNegative) -> std::tuple<std::vector<Point3D>, std::vector<TriangleMesh>, std::vector<std::pair<Point3D, Point3D>>, std::vector<bool>> {
        
        std::vector<Point3D> zeroVerts;
        std::vector<TriangleMesh> triangles;
        std::vector<std::pair<Point3D, Point3D>> zeroCrossingEdges; // Edges for each zero-crossing (as Point3D pairs)
        std::vector<bool> usedMidpoint; // Boolean mask: true if zero-crossing used midpoint interpolation
        

        zeroVerts.reserve(4*tets.size());
        triangles.reserve(2*tets.size());
        zeroCrossingEdges.reserve(4*tets.size());
        usedMidpoint.reserve(4*tets.size());


        // Build mapping from snappable vertices to their new indices
        std::vector<int> newInds(vts.size(), -1);
        int ct = 0;
        for (size_t i = 0; i < vts.size(); ++i) {
            if (mask[i]) {
                newInds[i] = ct++;
            } else {
                newInds[i] = ct;
            }
        }
        int numSnappable = ct;
        
        // Add snappable vertices first
        for (size_t i = 0; i < vts.size(); ++i) {
            if (mask[i]) {
                zeroVerts.push_back(vts[i]);
                // For snapped vertices, use the vertex itself as both endpoints (not a real edge)
                zeroCrossingEdges.push_back({vts[i], vts[i]});
                usedMidpoint.push_back(false); // Snapped vertices don't use midpoint
            }
        }
        
        // Hash map for zero-crossing points
        std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeHash;
        edgeHash.reserve(tets.size() * 3);
        
        int zeroVertCounter = numSnappable;
        
        // Process each tet
        for (const auto& tet : tets) {
            int vIndices[4] = {tet.v1, tet.v2, tet.v3, tet.v4};
            
            int signs[4];
            for (int i = 0; i < 4; ++i) {
                signs[i] = getSign(vals[vIndices[i]]);
            }
            
            int patternIdx = signs[0] * 8 + signs[1] * 4 + signs[2] * 2 + signs[3] * 1;
            const auto& entry = tetTableEinds[patternIdx];
            
            if (entry.edgeIndices.empty()) {
                continue;
            }
            
            int ori = getTetOrientation(
                vts[vIndices[0]], vts[vIndices[1]],
                vts[vIndices[2]], vts[vIndices[3]]
            );
            
            std::vector<int> localEdgeToZeroVert(6, -1);
            
            for (int edgeIdx1Based : entry.edgeIndices) {
                int edgeIdx = edgeIdx1Based - 1;
                const auto& edgePair = tetEdges[edgeIdx];
                int vIdx1 = vIndices[edgePair.first];
                int vIdx2 = vIndices[edgePair.second];
                
                int sign1 = getSign(vals[vIdx1]);
                int sign2 = getSign(vals[vIdx2]);
                if (sign1 == sign2) continue;
                
                EdgeKey edgeKey(vIdx1, vIdx2);
                
                int h = -1;
                // Check if we should snap
                if (sign1 == 0 && mask[vIdx1]) {
                    h = newInds[vIdx1];
                } else if (sign2 == 0 && mask[vIdx2]) {
                    h = newInds[vIdx2];
                } else {
                    // Use zero-crossing
                    auto it = edgeHash.find(edgeKey);
                    if (it != edgeHash.end()) {
                        h = it->second;
                    } else {
                        Point3D zeroPt;
                        bool usedMid = false;
                        // If either vertex is a new refined vertex (index >= originalVertexCount),
                        // use midpoint instead of linear interpolation
                        if (((vIdx1 >= originalVertexCount) || (vIdx2 >= originalVertexCount))||(vIdx1 <= changedFromPositiveToNegative.size()&&changedFromPositiveToNegative[vIdx1]) || (vIdx2 <= changedFromPositiveToNegative.size()&&changedFromPositiveToNegative[vIdx2])) {
                            // Use midpoint
                            usedMid = true;
                            zeroPt = Point3D(
                                (vts[vIdx1].x + vts[vIdx2].x) * 0.5,
                                (vts[vIdx1].y + vts[vIdx2].y) * 0.5,
                                (vts[vIdx1].z + vts[vIdx2].z) * 0.5
                            );
                        } else {
                            // Use linear interpolation
                            usedMid = true;
                            zeroPt = getLinearRoot(
                                vts[vIdx1], vts[vIdx2],
                                vals[vIdx1], vals[vIdx2]
                            );
                        }
                        h = zeroVertCounter++;
                        zeroVerts.push_back(zeroPt);
                        zeroCrossingEdges.push_back({vts[vIdx1], vts[vIdx2]}); // Store the edge endpoints as Point3D
                        usedMidpoint.push_back(usedMid); // Store midpoint usage
                        edgeHash[edgeKey] = h;
                    }
                }
                localEdgeToZeroVert[edgeIdx] = h;
            }
            
            // Create triangles
            for (const auto& cycle : entry.triangles) {
                std::vector<int> triVerts;
                for (int edgeIdx1Based : cycle) {
                    int edgeIdx = edgeIdx1Based - 1;
                    int zeroVertIdx = localEdgeToZeroVert[edgeIdx];
                    if (zeroVertIdx >= 0) {
                        triVerts.push_back(zeroVertIdx);
                    }
                }
                
                if (triVerts.size() == 3) {
                    // Check for duplicates
                    if (triVerts[0] != triVerts[1] && triVerts[1] != triVerts[2] && triVerts[0] != triVerts[2]) {
                        if (ori == 1) {
                            triangles.emplace_back(triVerts[0], triVerts[1], triVerts[2]);
                        } else {
                            triangles.emplace_back(triVerts[2], triVerts[1], triVerts[0]);
                        }
                    }
                }
            }
        }
        
        zeroVerts.shrink_to_fit();
        triangles.shrink_to_fit();
        zeroCrossingEdges.shrink_to_fit();
        usedMidpoint.shrink_to_fit();

        return {zeroVerts, triangles, zeroCrossingEdges, usedMidpoint};
    };
    
    auto refinementEndTime = std::chrono::steady_clock::now();
    auto refinementDuration = std::chrono::duration_cast<std::chrono::milliseconds>(refinementEndTime - refinementStartTime).count();
    std::cout << "Refinement took " << refinementDuration << " ms" << std::endl;
    auto pipelineDuration = std::chrono::duration_cast<std::chrono::milliseconds>(refinementEndTime - pipelineStartTime).count();
    std::cout << "  [PROFILE] refine_pipeline_total: " << pipelineDuration << " ms" << std::endl;

    // "isosurface time (refinement time)" timing for JSONL export.
    isosurfaceTimeMs = static_cast<double>(refinementDuration)/2.0;
    appendJsonlMetrics(simplex2TimeMs, isosurfaceTimeMs);
    auto marchTetSnapStartTime = std::chrono::steady_clock::now();
    auto [snapVerts, snapTris, snapZeroCrossingEdges, snapUsedMidpoint] = marchTetSnap(refinedVerts, refinedTetElements, refinedVals, snapMask, n, changedFromPositiveToNegative);
    std::cout << "Marching tetrahedra with snapping produced " << snapVerts.size()
              << " vertices and " << snapTris.size() << " triangles" << std::endl;

    //turn snapTris into a vector of vert indices
    std::vector<std::vector<int>> snapTrisVerts;
    for (const auto& t : snapTris) {
        snapTrisVerts.push_back({t.v1, t.v2, t.v3});
    }

    // Compute non-manifold edges and vertices
    nmeds = getSimplicialNonmaniEdges(snapTrisVerts);
    nmvts = getPolygonalNonmaniVerts(snapTrisVerts);
    
    std::cout << "Found " << nmeds.size() << " non-manifold exterior edges after snapping" << std::endl;
    std::cout << "Found " << nmvts.size() << " non-manifold exterior vertices after snapping" << std::endl;

    reportMeshManifoldness(snapVerts, snapTris, "Marching Tetrahedra with Snapping (after topoMinCut)");
    
    // Write result
    {
        std::filesystem::create_directories("output");
        std::ofstream ofs("output/new_surface_before_fairing.obj");
        ofs << std::fixed << std::setprecision(15);
        
        for (const auto& v : snapVerts) {
            ofs << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        
        for (const auto& t : snapTris) {
            ofs << "f " << (t.v1 + 1) << " " << (t.v2 + 1) << " " << (t.v3 + 1) << "\n";
        }
        
        ofs.close();
    }
    std::cout << "Wrote marching tetrahedra with snapping result to output/new_surface_before_fairing.obj" << std::endl;
    
    
    // Fairing algorithm
    {
        // Step 1: Build incidentVerts - for each vertex, record indices of vertices that share a triangle with it
        std::vector<std::vector<int>> incidentVerts(snapVerts.size());
        for (const auto& tri : snapTrisVerts) {
            // Each triangle has 3 vertices, so each vertex is incident to the other 2
            for (int i = 0; i < 3; ++i) {
                int vIdx = tri[i];
                for (int j = 0; j < 3; ++j) {
                    if (i != j) {
                        int otherVIdx = tri[j];
                        incidentVerts[vIdx].push_back(otherVIdx);
                    }
                }
            }
        }
        
        // Remove duplicates by converting to set and back to vector
        for (auto& vertList : incidentVerts) {
            std::unordered_set<int> uniqueVerts(vertList.begin(), vertList.end());
            vertList.assign(uniqueVerts.begin(), uniqueVerts.end());
        }
        
        // Step 2: Transform snapZeroCrossingEdges - pad towards center (10% on each end)
        std::vector<std::pair<Point3D, Point3D>> transformedEdges = snapZeroCrossingEdges;
        for (auto& edge : transformedEdges) {
            Point3D dir(
                edge.second.x - edge.first.x,
                edge.second.y - edge.first.y,
                edge.second.z - edge.first.z
            );
            // Shrink by 10% on each end
            Point3D offset(
                dir.x * 0.1,
                dir.y * 0.1,
                dir.z * 0.1
            );
            edge.first.x += offset.x;
            edge.first.y += offset.y;
            edge.first.z += offset.z;
            edge.second.x -= offset.x;
            edge.second.y -= offset.y;
            edge.second.z -= offset.z;
        }
        
        // Helper function: Project point onto line segment and return closest point
        auto projectOntoLineSegment = [](const Point3D& p, const Point3D& segStart, const Point3D& segEnd) -> Point3D {
            Point3D segDir(
                segEnd.x - segStart.x,
                segEnd.y - segStart.y,
                segEnd.z - segStart.z
            );
            
            Point3D toPoint(
                p.x - segStart.x,
                p.y - segStart.y,
                p.z - segStart.z
            );
            
            double segLenSq = segDir.x * segDir.x + segDir.y * segDir.y + segDir.z * segDir.z;
            if (segLenSq < 1e-12) {
                // Degenerate segment, return start point
                return segStart;
            }
            
            double t = (toPoint.x * segDir.x + toPoint.y * segDir.y + toPoint.z * segDir.z) / segLenSq;
            
            // Clamp to [0, 1] to stay on segment
            t = std::max(0.0, std::min(1.0, t));
            
            return Point3D(
                segStart.x + t * segDir.x,
                segStart.y + t * segDir.y,
                segStart.z + t * segDir.z
            );
        };
        
        // Step 3: Iterative fairing (20 iterations)
        const int numIterations = 20;
        for (int iter = 0; iter < numIterations; ++iter) {
            std::vector<Point3D> newVerts = snapVerts; // Temporary storage
            
            // Process each vertex that used midpoint
            for (size_t i = 0; i < snapVerts.size(); ++i) {
                if (snapUsedMidpoint[i]) {
                    // Compute 3D centroid of incident vertices
                    if (incidentVerts[i].empty()) {
                        continue; // Skip if no incident vertices
                    }
                    
                    Point3D centroid(0, 0, 0);
                    for (int j : incidentVerts[i]) {
                        centroid.x += snapVerts[j].x;
                        centroid.y += snapVerts[j].y;
                        centroid.z += snapVerts[j].z;
                    }
                    double count = static_cast<double>(incidentVerts[i].size());
                    centroid.x /= count;
                    centroid.y /= count;
                    centroid.z /= count;
                    
                    // Project centroid onto transformed zero-crossing edge
                    const auto& edge = transformedEdges[i];
                    Point3D projected = projectOntoLineSegment(centroid, edge.first, edge.second);
                    
                    // Store new coordinate
                    newVerts[i] = projected;
                }
            }
            
            // Update all vertices at once
            snapVerts = newVerts;
        }
        
        std::cout << "Fairing complete: " << numIterations << " iterations applied" << std::endl;
    }
    

    // Write result
    {
        std::filesystem::create_directories("output");
        std::ofstream ofs("output/new_surface_after_fairing.obj");
        ofs << std::fixed << std::setprecision(15);
        
        for (const auto& v : snapVerts) {
            ofs << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        
        for (const auto& t : snapTris) {
            ofs << "f " << (t.v1 + 1) << " " << (t.v2 + 1) << " " << (t.v3 + 1) << "\n";
        }
        
        ofs.close();
    }
    std::cout << "Wrote marching tetrahedra with snapping result to output/new_surface_after_fairing.obj" << std::endl;
    
    // ============================================================================
    // END REFINEMENT AND SNAPPING IMPLEMENTATION
    // ============================================================================

    
    // Option B cleanup: disable object-dependent debug exports on array-only path.
    std::cout << "Skipping legacy object-based WL/PLY exports on array-only path." << std::endl;
#if 0
    // Output all cells with updated alpha values to new_all_cells.wl
    std::cout << "Writing all cells with updated alpha values to Mathematica file..." << std::endl;
    std::filesystem::create_directories("output");
    {
        std::ofstream ofs("output/new_all_cells.wl");
        // Set fixed-point notation with sufficient precision to avoid scientific notation
        // Mathematica interprets 'e' as the constant E, so we avoid scientific notation
        ofs << std::fixed << std::setprecision(15);
        
        ofs << "(* All cells in order: {vertices, edges, triangles, tets} with coordinates and updated alpha values *)\n";
        ofs << "(* Format: {cellIndex, dimension, coordinates, alphaValue} *)\n";
        ofs << "(* Alpha values: updated from allSimplices after TopoMinCut optimization *)\n";
        ofs << "allCells = {\n";
        
        // Output vertices (dimension 0)
        ofs << "  (* Vertices (dimension 0) *)\n  {\n";
        bool firstVertex = true;
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            
            if (s->dimension == 0) {
                Vertex* v = static_cast<Vertex*>(s);
                if (!firstVertex) ofs << ",\n";
                ofs << "    {" << s->id << ", 0, {" << v->x << ", " << v->y << ", " << v->z 
                    << "}, " << s->val << "}";
                firstVertex = false;
            }
        }
        ofs << "\n  },\n";
        
        // Output edges (dimension 1)
        ofs << "  (* Edges (dimension 1) *)\n  {\n";
        bool firstEdge = true;
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            if (s->dimension == 1) {
                Edge* e = static_cast<Edge*>(s);
                if (!firstEdge) ofs << ",\n";
                ofs << "    {" << s->id << ", 1, {";
                for (size_t j = 0; j < e->verts.size(); ++j) {
                    if (j > 0) ofs << ", ";
                    ofs << "{" << e->verts[j][0] << ", " << e->verts[j][1] << ", " << e->verts[j][2] << "}";
                }
                ofs << "}, " << s->val << "}";
                firstEdge = false;
            }
        }
        ofs << "\n  },\n";
        
        // Output triangles (dimension 2)
        ofs << "  (* Triangles (dimension 2) *)\n  {\n";
        bool firstTriangle = true;
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            if (s->dimension == 2) {
                Triangle* t = static_cast<Triangle*>(s);
                if (!firstTriangle) ofs << ",\n";
                ofs << "    {" << s->id << ", 2, {";
                for (size_t j = 0; j < t->verts.size(); ++j) {
                    if (j > 0) ofs << ", ";
                    ofs << "{" << t->verts[j][0] << ", " << t->verts[j][1] << ", " << t->verts[j][2] << "}";
                }
                ofs << "}, " << s->val << "}";
                firstTriangle = false;
            }
        }
        ofs << "\n  },\n";
        
        // Output tets (dimension 3)
        ofs << "  (* Tetrahedra (dimension 3) *)\n  {\n";
        bool firstTet = true;
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            if (s->dimension == 3) {
                Tet* tet = static_cast<Tet*>(s);
                if (!firstTet) ofs << ",\n";
                ofs << "    {" << s->id << ", 3, {";
                for (size_t j = 0; j < tet->verts.size(); ++j) {
                    if (j > 0) ofs << ", ";
                    ofs << "{" << tet->verts[j][0] << ", " << tet->verts[j][1] << ", " << tet->verts[j][2] << "}";
                }
                ofs << "}, " << s->val << "}";
                firstTet = false;
            }
        }
        ofs << "\n  }\n";
        
        ofs << "};\n\n";
        ofs << "(* Helper functions for visualization *)\n";
        ofs << "getVertices[] := allCells[[1]];\n";
        ofs << "getEdges[] := allCells[[2]];\n";
        ofs << "getTriangles[] := allCells[[3]];\n";
        ofs << "getTets[] := allCells[[4]];\n\n";
        ofs << "(* Example Manipulate visualization with alpha filtering *)\n";
        ofs << "Manipulate[\n";
        ofs << "  Module[{filteredVertices, filteredEdges, filteredTriangles, filteredTets},\n";
        ofs << "    (* Filter cells by alpha value *)\n";
        ofs << "    filteredVertices = Select[getVertices[], alphaMin <= #[[4]] <= alphaMax &];\n";
        ofs << "    filteredEdges = Select[getEdges[], alphaMin <= #[[4]] <= alphaMax &];\n";
        ofs << "    filteredTriangles = Select[getTriangles[], alphaMin <= #[[4]] <= alphaMax &];\n";
        ofs << "    filteredTets = Select[getTets[], alphaMin <= #[[4]] <= alphaMax &];\n";
        ofs << "    \n";
        ofs << "    Graphics3D[{\n";
        ofs << "      (* Vertices *)\n";
        ofs << "      If[showVertices && Length[filteredVertices] > 0, {PointSize[0.01], Red, Point[filteredVertices[[All, 3]]]}, {}],\n";
        ofs << "      (* Edges *)\n";
        ofs << "      If[showEdges && Length[filteredEdges] > 0, {Thickness[0.001], Blue, Line[filteredEdges[[All, 3]]]}, {}],\n";
        ofs << "      (* Triangles *)\n";
        ofs << "      If[showTriangles && Length[filteredTriangles] > 0, {Opacity[0.05], Green, EdgeForm[], Polygon[filteredTriangles[[All, 3]]]}, {}],\n";
        ofs << "      (* Tets - shown as faces *)\n";
        ofs << "      If[showTets && Length[filteredTets] > 0, {Opacity[0.2], Yellow, EdgeForm[{Thin, Black}], \n";
        ofs << "        Flatten[Table[Polygon[Subsets[filteredTets[[i, 3]], {3}]], {i, Length[filteredTets]}]]}, {}]\n";
        ofs << "    }, Boxed -> True, Axes -> True, ImageSize -> Large, PlotLabel -> \n";
        ofs << "      StringForm[\"Showing: `` vertices, `` edges, `` triangles, `` tets\", \n";
        ofs << "        Length[filteredVertices], Length[filteredEdges], Length[filteredTriangles], Length[filteredTets]]]\n";
        ofs << "  ],\n";
        ofs << "  {{showVertices, False, \"Show Vertices\"}, {True, False}},\n";
        ofs << "  {{showEdges, True, \"Show Edges\"}, {True, False}},\n";
        ofs << "  {{showTriangles, True, \"Show Triangles\"}, {True, False}},\n";
        ofs << "  {{showTets, False, \"Show Tets\"}, {True, False}},\n";
        ofs << "  {{alphaMin, Min[Flatten[allCells[[All, All, 4]]]], \"Alpha Min\"}, Min[Flatten[allCells[[All, All, 4]]]], Max[Flatten[allCells[[All, All, 4]]]]},\n";
        ofs << "  {{alphaMax, Max[Flatten[allCells[[All, All, 4]]]], \"Alpha Max\"}, Min[Flatten[allCells[[All, All, 4]]]], Max[Flatten[allCells[[All, All, 4]]]]}\n";
        ofs << "]\n";
        ofs.close();
    }
#endif
    // Legacy object-dependent exports intentionally disabled on array-only path.
#if 0
    std::cout << "Wrote all cells with updated alpha values to output/new_all_cells.wl" << std::endl;

    // Output faces with 2 incident tets where exactly one tet has alpha <= 0 from allSimplices
    {
        std::ofstream ofs("output/new_boundary_faces_alpha_leq_0.wl");
        ofs << std::fixed << std::setprecision(15);
        ofs << "(* Faces with 2 incident tets where exactly one tet has alpha <= 0 *)\n";
        ofs << "(* Format: Red polygons *)\n";
        ofs << "Graphics3D[{\n";
        ofs << "  Red, Opacity[0.5],\n";
        ofs << "  Polygon[{\n";
        bool first = true;
        
        // Build mapping from triangle to incident tets (single pass through tets)
        std::unordered_map<Triangle*, std::vector<Tet*>> triangleToTets;
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            if (s->dimension == 3) {
                Tet* tet = static_cast<Tet*>(s);
                for (Simplex* face : tet->faces) {
                    if (face->dimension == 2) {
                        Triangle* t = static_cast<Triangle*>(face);
                        triangleToTets[t].push_back(tet);
                    }
                }
            }
        }
        
        // Now iterate through triangles and check the mapping
        for (size_t i = 0; i < allSimplices.size(); ++i) {
            Simplex* s = allSimplices[i];
            if (s->dimension == 2) {
                Triangle* t = static_cast<Triangle*>(s);
                auto it = triangleToTets.find(t);
                if (it != triangleToTets.end() && it->second.size() == 2) {
                    bool tet1_leq0 = it->second[0]->val <= 0.0;
                    bool tet2_leq0 = it->second[1]->val <= 0.0;
                    // Exactly one tet has alpha <= 0
                    if (tet1_leq0 != tet2_leq0) {
                        if (!first) ofs << ",\n";
                        ofs << "    {";
                        for (size_t j = 0; j < t->verts.size(); ++j) {
                            if (j > 0) ofs << ", ";
                            ofs << "{" << t->verts[j][0] << ", " << t->verts[j][1] << ", " << t->verts[j][2] << "}";
                        }
                        ofs << "}";
                        first = false;
                    }
                }
            }
        }
        ofs << "\n  }]\n";
        ofs << "}, Boxed -> True, Axes -> True, ImageSize -> Large]\n";
        ofs.close();
    }
    std::cout << "Wrote boundary faces with alpha <= 0 to output/new_boundary_faces_alpha_leq_0.wl" << std::endl;

    //instead of reading a file, I want to read all files in a directory and store them in a vector of vectors
    //so do a for loop

    
    std::vector<int> remaining_negative = tet_topo_out.remaining_negative;
    std::vector<int> protected_indices = tet_topo_out.protected_indices;

    // Output remaining_negative to PLY format
    {
        std::cout << "Writing remaining_negative to PLY files..." << std::endl;
        
        // Ensure output directory exists
        std::filesystem::create_directories("output");
        
        // First pass: collect all unique vertices
        std::map<std::tuple<double, double, double>, int> vertexMap;
        std::vector<std::tuple<double, double, double>> vertices;
        std::vector<std::pair<int, int>> edges;
        std::vector<std::vector<int>> faces; // Store face vertex indices
        int vertexIndex = 0; // PLY uses 0-based indexing
        
        // Collect vertices from dimension 0 simplices
        for (int id : remaining_negative) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 0) {
                    std::tuple<double, double, double> vertex = std::make_tuple(s->x, s->y, s->z);
                    
                    if (vertexMap.find(vertex) == vertexMap.end()) {
                        vertexMap[vertex] = vertexIndex;
                        vertices.push_back(vertex);
                        vertexIndex++;
                    }
            }
        }
        
        // Collect vertices from edges (dimension 1)
        for (int id : remaining_negative) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 1) {
                    Edge* edge = static_cast<Edge*>(s);
                    if (edge->verts.size() >= 2) {
                        std::tuple<double, double, double> v1 = std::make_tuple(
                            edge->verts[0][0], edge->verts[0][1], edge->verts[0][2]);
                        std::tuple<double, double, double> v2 = std::make_tuple(
                            edge->verts[1][0], edge->verts[1][1], edge->verts[1][2]);
                        
                        if (vertexMap.find(v1) == vertexMap.end()) {
                            vertexMap[v1] = vertexIndex;
                            vertices.push_back(v1);
                            vertexIndex++;
                        }
                        if (vertexMap.find(v2) == vertexMap.end()) {
                            vertexMap[v2] = vertexIndex;
                            vertices.push_back(v2);
                            vertexIndex++;
                        }
                        // Create edge
                        edges.push_back(std::make_pair(vertexMap[v1], vertexMap[v2]));
                    }
            }
        }
        
        // Process triangles and collect any missing vertices, then create faces
        for (int id : remaining_negative) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 2) {
                    Triangle* triangle = static_cast<Triangle*>(s);
                    if (triangle->verts.size() >= 3) {
                        std::vector<std::tuple<double, double, double>> triVertices;
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[0][0], triangle->verts[0][1], triangle->verts[0][2]));
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[1][0], triangle->verts[1][1], triangle->verts[1][2]));
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[2][0], triangle->verts[2][1], triangle->verts[2][2]));
                        
                        // Add any missing vertices and get indices
                        std::vector<int> indices;
                        for (const auto& v : triVertices) {
                            if (vertexMap.find(v) == vertexMap.end()) {
                                vertexMap[v] = vertexIndex;
                                vertices.push_back(v);
                                vertexIndex++;
                            }
                            indices.push_back(vertexMap[v]);
                        }
                        
                        // Store triangle as a face with 3 indices
                        if (indices.size() == 3) {
                            faces.push_back(indices);
                        }
                    }
            }
        }
        
        // Write PLY file: edges only (ASCII format, 0-based indexing)
        std::ofstream plyEdgesFile("output/remaining_negative_edgesOnly.ply");
        plyEdgesFile << "ply\n";
        plyEdgesFile << "format ascii 1.0\n";
        plyEdgesFile << "comment TopoTet output\n";
        plyEdgesFile << "element vertex " << vertices.size() << "\n";
        plyEdgesFile << "property float x\n";
        plyEdgesFile << "property float y\n";
        plyEdgesFile << "property float z\n";
        plyEdgesFile << "element edge " << edges.size() << "\n";
        plyEdgesFile << "property int vertex1\n";
        plyEdgesFile << "property int vertex2\n";
        plyEdgesFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyEdgesFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }
        
        // Write edges (0-based)
        for (const auto& e : edges) {
            plyEdgesFile << e.first << " " << e.second << "\n";
        }
        
        plyEdgesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << edges.size() 
                  << " edges to output/remaining_negative_edgesOnly.ply" << std::endl;
        
        // Write PLY file: faces only (ASCII format, 0-based indexing)
        std::ofstream plyFacesFile("output/remaining_negative_facesOnly.ply");
        plyFacesFile << "ply\n";
        plyFacesFile << "format ascii 1.0\n";
        plyFacesFile << "comment TopoTet output\n";
        plyFacesFile << "element vertex " << vertices.size() << "\n";
        plyFacesFile << "property float x\n";
        plyFacesFile << "property float y\n";
        plyFacesFile << "property float z\n";
        plyFacesFile << "element face " << faces.size() << "\n";
        plyFacesFile << "property list uchar int vertex_indices\n";
        plyFacesFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyFacesFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }
        
        // Write faces (0-based)
        for (const auto& face : faces) {
            plyFacesFile << face.size();
            for (int idx : face) {
                plyFacesFile << " " << idx;
            }
            plyFacesFile << "\n";
        }
        plyFacesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << faces.size() 
                  << " faces to output/remaining_negative_facesOnly.ply" << std::endl;
    }
    
    // Output protected_indices to PLY format
    {
        std::cout << "Writing protected_indices to PLY files..." << std::endl;
        
        // Ensure output directory exists
        std::filesystem::create_directories("output");
        
        // First pass: collect all unique vertices
        std::map<std::tuple<double, double, double>, int> vertexMap;
        std::vector<std::tuple<double, double, double>> vertices;
        std::vector<std::pair<int, int>> edges;
        std::vector<std::vector<int>> faces; // Store face vertex indices
        int vertexIndex = 0; // PLY uses 0-based indexing
        
        // Collect vertices from dimension 0 simplices
        for (int id : protected_indices) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 0) {
                    std::tuple<double, double, double> vertex = std::make_tuple(s->x, s->y, s->z);
                    
                    if (vertexMap.find(vertex) == vertexMap.end()) {
                        vertexMap[vertex] = vertexIndex;
                        vertices.push_back(vertex);
                        vertexIndex++;
                    }
            }
        }
        
        // Collect vertices from edges (dimension 1)
        for (int id : protected_indices) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 1) {
                    Edge* edge = static_cast<Edge*>(s);
                    if (edge->verts.size() >= 2) {
                        std::tuple<double, double, double> v1 = std::make_tuple(
                            edge->verts[0][0], edge->verts[0][1], edge->verts[0][2]);
                        std::tuple<double, double, double> v2 = std::make_tuple(
                            edge->verts[1][0], edge->verts[1][1], edge->verts[1][2]);
                        
                        if (vertexMap.find(v1) == vertexMap.end()) {
                            vertexMap[v1] = vertexIndex;
                            vertices.push_back(v1);
                            vertexIndex++;
                        }
                        if (vertexMap.find(v2) == vertexMap.end()) {
                            vertexMap[v2] = vertexIndex;
                            vertices.push_back(v2);
                            vertexIndex++;
                        }
                        // Create edge
                        edges.push_back(std::make_pair(vertexMap[v1], vertexMap[v2]));
                    }
            }
        }
        
        // Process triangles and collect any missing vertices, then create faces
        for (int id : protected_indices) {
            if (Simplex* s = getKeptSimplex(id); s != nullptr && s->dimension == 2) {
                    Triangle* triangle = static_cast<Triangle*>(s);
                    if (triangle->verts.size() >= 3) {
                        std::vector<std::tuple<double, double, double>> triVertices;
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[0][0], triangle->verts[0][1], triangle->verts[0][2]));
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[1][0], triangle->verts[1][1], triangle->verts[1][2]));
                        triVertices.push_back(std::make_tuple(
                            triangle->verts[2][0], triangle->verts[2][1], triangle->verts[2][2]));
                        
                        // Add any missing vertices and get indices
                        std::vector<int> indices;
                        for (const auto& v : triVertices) {
                            if (vertexMap.find(v) == vertexMap.end()) {
                                vertexMap[v] = vertexIndex;
                                vertices.push_back(v);
                                vertexIndex++;
                            }
                            indices.push_back(vertexMap[v]);
                        }
                        
                        // Store triangle as a face with 3 indices
                        if (indices.size() == 3) {
                            faces.push_back(indices);
                        }
                    }
            }
        }
        
        // Write PLY file: edges only (ASCII format, 0-based indexing)
        std::ofstream plyEdgesFile("output/original_skeleton_edgesOnly.ply");
        plyEdgesFile << "ply\n";
        plyEdgesFile << "format ascii 1.0\n";
        plyEdgesFile << "comment TopoTet output\n";
        plyEdgesFile << "element vertex " << vertices.size() << "\n";
        plyEdgesFile << "property float x\n";
        plyEdgesFile << "property float y\n";
        plyEdgesFile << "property float z\n";
        plyEdgesFile << "element edge " << edges.size() << "\n";
        plyEdgesFile << "property int vertex1\n";
        plyEdgesFile << "property int vertex2\n";
        plyEdgesFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyEdgesFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }
        
        // Write edges (0-based)
        for (const auto& e : edges) {
            plyEdgesFile << e.first << " " << e.second << "\n";
        }
        
        plyEdgesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << edges.size() 
                  << " edges to output/original_skeleton_edgesOnly.ply" << std::endl;
        
        // Write PLY file: faces only (ASCII format, 0-based indexing)
        std::ofstream plyFacesFile("output/original_skeleton_facesOnly.ply");
        plyFacesFile << "ply\n";
        plyFacesFile << "format ascii 1.0\n";
        plyFacesFile << "comment TopoTet output\n";
        plyFacesFile << "element vertex " << vertices.size() << "\n";
        plyFacesFile << "property float x\n";
        plyFacesFile << "property float y\n";
        plyFacesFile << "property float z\n";
        plyFacesFile << "element face " << faces.size() << "\n";
        plyFacesFile << "property list uchar int vertex_indices\n";
        plyFacesFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyFacesFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }
        
        // Write faces (0-based)
        for (const auto& face : faces) {
            plyFacesFile << face.size();
            for (int idx : face) {
                plyFacesFile << " " << idx;
            }
            plyFacesFile << "\n";
        }
        plyFacesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << faces.size() 
                  << " faces to output/original_skeleton_facesOnly.ply" << std::endl;
    }
#endif

    // Clean up memory
    for (Edge* edge : edges) {
        delete edge;
    }
    for (Triangle* triangle : triangles) {
        delete triangle;
    }

    std::cout << "Simplice generation completed in " << simplice_generation_duration.count() << " milliseconds" << std::endl;
    
    std::cout << "Processing completed successfully" << std::endl;
    return 0;

}

} // namespace tet
} // namespace prescribed_topo

