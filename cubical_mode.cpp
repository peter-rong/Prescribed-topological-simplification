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
#include <filesystem>
#include <unordered_set>
#include <map>
#ifdef _WIN32
#include <io.h>
#endif
#include "mrc_reader.h"
#include "marching_cubes.h"

#include <Eigen/Sparse>
#include "topomincut_runner.hpp"

namespace prescribed_topo {
namespace cubical {

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
                                   const std::vector<Triangle>& triangles,
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
class Quad;
class Cube;


// Base class for all simplices (memory optimized)
class Simplex {
public:
    int id;
    double val;
    int dimension;
    int x,y,z;
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
    
    virtual std::vector<Simplex*> getFaces() const {
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
    
    Vertex(int id, double val, int x, int y, int z) {
        this->id = id;
        this->val = val;
        this->dimension = 0;
        this->x = x;
        this->y = y;
        this->z = z;
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
    std::vector<Quad*> incident_quads;
    
    Edge(int id, double val, int x, int y, int z, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 1;
        this->faces = children;
        this->x = x;
        this->y = y;
        this->z = z;

        // Pre-allocate incident_quads
        incident_quads.reserve(4); // Max 4 quads per edge in 3D grid
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
        incident_quads.shrink_to_fit();
    }
};

class Quad : public Simplex {
public:
    std::vector<Cube*> incident_cubes;
    
    Quad(int id, double val, int x, int y, int z, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 2;
        this->faces = children;
        this->x = x;
        this->y = y;
        this->z = z;
        // Pre-allocate incident_cubes
        incident_cubes.reserve(2); // Max 2 cubes per quad
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
        incident_cubes.shrink_to_fit();
    }
};

class Cube : public Simplex {
public:
    
    Cube(int id, double val, int x, int y, int z, const std::vector<Simplex*>& children) {
        this->id = id;
        this->val = val;
        this->dimension = 3;
        this->faces = children;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    void optimizeMemory() override {
        Simplex::optimizeMemory();
    }
};

// Filter result for keeping only in-range ∪ interface simplices (memory optimized)
struct FilterResult {
    std::vector<Simplex*> keptSimplices;               // kept simplices in new order
    std::vector<std::vector<int>> keptChildrenByCol;   // child rows per kept column (0-based), new indices
    std::vector<int> oldToNew;                         // size N, -1 if dropped
    std::vector<int> newToOld;     
                        // size K
    void deleteAll() {
        keptSimplices.clear();
        keptSimplices.shrink_to_fit();
        for (auto& col : keptChildrenByCol) {
            col.clear();
            col.shrink_to_fit();
        }
        keptChildrenByCol.clear();
        keptChildrenByCol.shrink_to_fit();
        oldToNew.clear();
        oldToNew.shrink_to_fit();
        newToOld.clear();
        newToOld.shrink_to_fit();
    }
    // Memory optimization: shrink vectors after construction
    void optimizeMemory() {
        keptSimplices.shrink_to_fit();
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
    const std::vector<Simplex*>& allSimplices,
    const std::vector<std::vector<int>>& childrenByColumn,
    double alpha1,
    double alpha2
    )
{
    const int n = static_cast<int>(allSimplices.size());
    std::vector<bool> inRange(n, false); // Use bool instead of char for memory efficiency
    for (int i = 0; i < n; ++i) {
        const double a = allSimplices[i]->val;
        
        inRange[i] = true;
    }

    //print how many are inRange for each dimension
    int totalInRangeVertices = 0;
    int totalInRangeEdges = 0;
    int totalInRangeQuads = 0;
    int totalInRangeCubes = 0;
    for (int i = 0; i < n; ++i) {
        if (inRange[i]){
            if (allSimplices[i]->dimension == 0) totalInRangeVertices++;
            if (allSimplices[i]->dimension == 1) totalInRangeEdges++;
            if (allSimplices[i]->dimension == 2) totalInRangeQuads++;
            if (allSimplices[i]->dimension == 3) totalInRangeCubes++;
        }
    }
    std::cout << "Total inRange vertices: " << totalInRangeVertices << std::endl;
    std::cout << "Total inRange edges: " << totalInRangeEdges << std::endl;
    std::cout << "Total inRange quads: " << totalInRangeQuads << std::endl;
    std::cout << "Total inRange cubes: " << totalInRangeCubes << std::endl;

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
    std::cout << "HERE" << std::endl;

    std::cout << "First index: " << firstIndex << std::endl;
    std::cout << "Last index: " << lastIndex << std::endl;

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

    parentsByRow.clear();
    parentsByRow.shrink_to_fit();

    //output how many are interface
    int totalInterface = 0;
    int totalInterfaceVertices = 0;
    int totalInterfaceEdges = 0;
    int totalInterfaceQuads = 0;
    int totalInterfaceCubes = 0;

    for (int i = 0; i < isInterface.size(); ++i) {
        if (isInterface[i]){ 
            
            totalInterface++;

            if (allSimplices[i]->dimension == 0) totalInterfaceVertices++;
            if (allSimplices[i]->dimension == 1) totalInterfaceEdges++;
            if (allSimplices[i]->dimension == 2) totalInterfaceQuads++;
            if (allSimplices[i]->dimension == 3) totalInterfaceCubes++;
        }
    }

    for (int i = lastIndex; i < n; ++i) {
        if (inRange[i]) continue;
        for (int c : childrenByColumn[i]) { if (inRange[c]||isInterface[c]) { isInterface[i] = true; break; } }

    }

    totalInterface = 0;
    totalInterfaceVertices = 0;
    totalInterfaceEdges = 0;
    totalInterfaceQuads = 0;
    totalInterfaceCubes = 0;

    for (int i = 0; i < isInterface.size(); ++i) {
        if (isInterface[i]){ 
            
            totalInterface++;

            if (allSimplices[i]->dimension == 0) totalInterfaceVertices++;
            if (allSimplices[i]->dimension == 1) totalInterfaceEdges++;
            if (allSimplices[i]->dimension == 2) totalInterfaceQuads++;
            if (allSimplices[i]->dimension == 3) totalInterfaceCubes++;
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

    // Kept simplices
    std::vector<Simplex*> keptSimplices;
    keptSimplices.reserve(k);
    for (int idx : newToOld) keptSimplices.push_back(allSimplices[idx]);

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

    FilterResult result{ std::move(keptSimplices), std::move(keptChildrenByCol), std::move(oldToNew), std::move(newToOld) };
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

// Structure to hold cell complex data
struct CellComplex {
    // cc[d][i] = list of indices of (d-1)-cells that bound the d-cell i
    // d = 0: vertices (no lower-dimensional faces, entries left empty)
    // d = 1: edges -> 2 vertex indices
    // d = 2: quads -> 4 edge indices
    // d = 3: cubes -> 6 quad indices
    std::vector<std::vector<std::vector<int>>> cc;       // size 4
    std::vector<std::vector<double>> vals;               // vals[d][i] = scalar value on that d-cell
    std::vector<Point3D> vertexCoords; // high-res grid coordinates for 0-cells
    
    CellComplex() : cc(4), vals(4) {}
};


// Function to convert high-res array to cell complex structure
std::pair<CellComplex, std::vector<bool>> createCellComplexFromHighRes(
    const std::vector<std::vector<std::vector<double>>>& highResArray,
    const std::vector<std::vector<std::vector<double>>>& toHalfVerts) {

    CellComplex result;
    std::vector<bool>vertNeedMidPoint;

    // Create padded highResArray: add a border of +100 in every direction
    const int hx = static_cast<int>(highResArray.size());
    const int hy = static_cast<int>(highResArray[0].size());
    const int hz = static_cast<int>(highResArray[0][0].size());

    // Single index array: same dimensions as highResArray, store index at each position (-1 if no cell)
    std::vector<std::vector<std::vector<int>>> indexArray(hx);
    
    // Initialize all to -1 (invalid index)
    for (int x = 0; x < hx; ++x) {
        indexArray[x].resize(hy);
        for (int y = 0; y < hy; ++y) {
            indexArray[x][y].resize(hz, -1);
        }
    }

    result.vals[0].reserve(hx * hy * hz);
    result.vals[1].reserve(hx * hy * hz);
    result.vals[2].reserve(hx * hy * hz);
    result.vals[3].reserve(hx * hy * hz);
    vertNeedMidPoint.reserve(hx * hy * hz);

    result.vertexCoords.reserve(hx * hy * hz);
    result.cc[0].reserve(hx * hy * hz);
    result.cc[1].reserve(hx * hy * hz);
    result.cc[2].reserve(hx * hy * hz);
    result.cc[3].reserve(hx * hy * hz);
    
    int vertexIdx = 0;
    int edgeIdx = 0;
    int quadIdx = 0;
    int cubeIdx = 0;

    for (int x = 0; x < hx-2; x += 2) {
        for (int y = 0; y < hy-2; y += 2) {
            for (int z = 0; z < hz-2; z += 2) {
            // Get the 8 vertex corners for this cube
                std::vector<std::tuple<int, int, int>> corners = {
                    {x,     y,     z},
                    {x+2,   y,     z},
                    {x,     y+2,   z},
                    {x+2,   y+2,   z},
                    {x,     y,     z+2},
                    {x+2,   y,     z+2},
                    {x,     y+2,   z+2},
                    {x+2,   y+2,   z+2}
                };

                std::vector<std::tuple<int, int, int>> edges = {
                    {x,   y+1,     z},
                    {x,     y,   z+1},
                    {x,     y+2,   z+1},
                    {x, y+1, z+2},
                    {x+2,   y+1,   z},
                    {x+2,   y,    z+1},
                    {x+2,   y+2,   z+1},
                    {x+2,   y+1,   z+2},
                    {x+1,     y,     z},
                    {x+1, y, z+2},
                    {x+1, y+2, z},
                    {x+1, y+2, z+2}

                };

                std::vector<std::tuple<int, int, int>> faces = {
                    {x+1,y+1, z},
                    {x+1,y+1, z+2},
                    {x+1, y, z+1},
                    {x+1, y+2, z+1},
                    {x, y+1, z+1},
                    {x+2,y+1, z+1}

                };

                std::vector<int> cube = {x+1,y+1,z+1};
                /*
                bool allPositive = true;
                for(const auto& [cx, cy, cz] : corners) {
                    // Defensive: Ensure inside bounds as this triple-nested loop is already restricted
                    if (highResArray[cx][cy][cz] <= 0) {
                        allPositive = false;
                        break;
                    }
                }
                if (allPositive) {
                    continue;
                }
                */
                // Check if all 27 values in the 3x3x3 cube from (x,y,z) to (x+2,y+2,z+2) are <= 0
                bool allNonPositive = true;
                for (int dx = 0; dx <= 2; ++dx) {
                    for (int dy = 0; dy <= 2; ++dy) {
                        for (int dz = 0; dz <= 2; ++dz) {
                            int nx = x + dx;
                            int ny = y + dy;
                            int nz = z + dz;
                            if (highResArray[nx][ny][nz] > 0) {
                                allNonPositive = false;
                                break;
                            }
                        }
                        if (!allNonPositive) break;
                    }
                    if (!allNonPositive) break;
                }

                bool allNonNegative= true;
                for (int dx = 0; dx <= 2; ++dx) {
                    for (int dy = 0; dy <= 2; ++dy) {
                        for (int dz = 0; dz <= 2; ++dz) {
                            int nx = x + dx;
                            int ny = y + dy;
                            int nz = z + dz;
                            if (highResArray[nx][ny][nz] <= 0) {
                                allNonNegative = false;
                                break;
                            }
                        }
                        if (!allNonNegative) break;
                    }
                    if (!allNonNegative) break;
                }
                
                if (allNonPositive || allNonNegative) {
                    continue;
                }
                
                //create vertices
                for(const auto& [cx, cy, cz] : corners) {

                    if (indexArray[cx][cy][cz] != -1) {
                        continue;
                    }
                    double val = highResArray[cx][cy][cz];
                    if (cx>=2 && cy>=2 && cz>=2 && cx<hx-2 && cy<hy-2 && cz<hz-2 && toHalfVerts[(cx-2)/2][(cy-2)/2][(cz-2)/2] == 1 ) {
                        vertNeedMidPoint.push_back(true);
                    }
                    else{
                        vertNeedMidPoint.push_back(false);
                    }   
                    // Vertex
                    indexArray[cx][cy][cz] = vertexIdx;
                    result.cc[0].push_back({vertexIdx});
                    vertexIdx += 1;
                    result.vertexCoords.push_back(Point3D(int(cx/2)-1, int(cy/2)-1, int(cz/2)-1));
                    result.vals[0].push_back(val);
                    
                }
                
                //create edges
                for(const auto& [cx, cy, cz] : edges) {
                    if (indexArray[cx][cy][cz] != -1) {
                        continue;
                    }
                    double val = highResArray[cx][cy][cz];
                    indexArray[cx][cy][cz] = edgeIdx;
                    
                    result.vals[1].push_back(val);

                    int xMod = cx % 2;
                    int yMod = cy % 2;
                    int zMod = cz % 2;
                    
                    std::vector<int> vertexIndices;
                    
                    if (xMod == 1) {

                        // x-edge: connects vertices at (x-1,y,z) and (x+1,y,z)
                        if (cx-1 >= 0 && indexArray[cx-1][cy][cz] != -1) {
                            vertexIndices.push_back(indexArray[cx-1][cy][cz]);
                        }
                        if (cx+1 < hx && indexArray[cx+1][cy][cz] != -1) {
                            vertexIndices.push_back(indexArray[cx+1][cy][cz]);
                        }
                    } else if (yMod == 1) {
                    
                        // y-edge: connects vertices at (x,y-1,z) and (x,y+1,z)
                        if (cy-1 >= 0 && indexArray[cx][cy-1][cz] != -1) {
                            vertexIndices.push_back(indexArray[cx][cy-1][cz]);
                        }
                        if (cy+1 < hy && indexArray[cx][cy+1][cz] != -1) {
                            vertexIndices.push_back(indexArray[cx][cy+1][cz]);
                        }
                    } else if (zMod == 1) {
                        
                        // z-edge: connects vertices at (x,y,z-1) and (x,y,z+1)
                        if (cz-1 >= 0 && indexArray[cx][cy][cz-1] != -1) {
                            vertexIndices.push_back(indexArray[cx][cy][cz-1]);
                        }
                        if (cz+1 < hz && indexArray[cx][cy][cz+1] != -1) {
                            vertexIndices.push_back(indexArray[cx][cy][cz+1]);
                        }
                    }
                    if (vertexIndices.size() != 2) {
                        std::cout << "edge has " << vertexIndices.size() << " vertices (error)" << std::endl;
                        continue;
                    }
                    result.cc[1].push_back(vertexIndices);
                    edgeIdx += 1;
                }

                
                //create faces
                for(const auto& [cx, cy, cz] : faces) {
                    if (indexArray[cx][cy][cz] != -1) {
                        continue;
                    }
                    double val = highResArray[cx][cy][cz];
                    indexArray[cx][cy][cz] = quadIdx;
                    
                    result.vals[2].push_back(val);

                    int xMod = cx % 2;
                    int yMod = cy % 2;
                    int zMod = cz % 2;

                    std::vector<int> edgeIndices;
                
                    if (xMod == 1 && yMod == 1) {
                        // xy-quad (perpendicular to z-axis): has 2 x-edges and 2 y-edges
                        if (cy-1 >= 0 && indexArray[cx][cy-1][cz] != -1) edgeIndices.push_back(indexArray[cx][cy-1][cz]);
                        if (cy+1 < hy && indexArray[cx][cy+1][cz] != -1) edgeIndices.push_back(indexArray[cx][cy+1][cz]);
                        if (cx-1 >= 0 && indexArray[cx-1][cy][cz] != -1) edgeIndices.push_back(indexArray[cx-1][cy][cz]);
                        if (cx+1 < hx && indexArray[cx+1][cy][cz] != -1) edgeIndices.push_back(indexArray[cx+1][cy][cz]);
                    } else if (xMod == 1 && zMod == 1) {
                        // xz-quad (perpendicular to y-axis): has 2 x-edges and 2 z-edges
                        if (cz-1 >= 0 && indexArray[cx][cy][cz-1] != -1) edgeIndices.push_back(indexArray[cx][cy][cz-1]);
                        if (cz+1 < hz && indexArray[cx][cy][cz+1] != -1) edgeIndices.push_back(indexArray[cx][cy][cz+1]);
                        if (cx-1 >= 0 && indexArray[cx-1][cy][cz] != -1) edgeIndices.push_back(indexArray[cx-1][cy][cz]);
                        if (cx+1 < hx && indexArray[cx+1][cy][cz] != -1) edgeIndices.push_back(indexArray[cx+1][cy][cz]);
                    } else if (yMod == 1 && zMod == 1) {
                        // yz-quad (perpendicular to x-axis): has 2 y-edges and 2 z-edges
                        if (cz-1 >= 0 && indexArray[cx][cy][cz-1] != -1) edgeIndices.push_back(indexArray[cx][cy][cz-1]);
                        if (cz+1 < hz && indexArray[cx][cy][cz+1] != -1) edgeIndices.push_back(indexArray[cx][cy][cz+1]);
                        if (cy-1 >= 0 && indexArray[cx][cy-1][cz] != -1) edgeIndices.push_back(indexArray[cx][cy-1][cz]);
                        if (cy+1 < hy && indexArray[cx][cy+1][cz] != -1) edgeIndices.push_back(indexArray[cx][cy+1][cz]);
                    }

                    if (edgeIndices.size() != 4) {
                        std::cout << "quad has " << edgeIndices.size() << " edges (error)" << std::endl;
                        continue;
                    }
                    
                    result.cc[2].push_back(edgeIndices);
                    quadIdx += 1;
                }

                int cx = cube[0];
                int cy = cube[1];
                int cz = cube[2];
                //create cube
                if (indexArray[cx][cy][cz] != -1) {
                    std::cout << "creating cube but already exists" << std::endl;
                    continue;
                }
                double val = highResArray[cx][cy][cz];
                indexArray[cx][cy][cz] = cubeIdx;

                std::vector<int> quadIndices;
                
                // Cube has 6 quads: one in each direction (x, y, z) at both ends
                if (cx-1 >= 0 && indexArray[cx-1][cy][cz] != -1) quadIndices.push_back(indexArray[cx-1][cy][cz]);
                if (cx+1 < hx && indexArray[cx+1][cy][cz] != -1) quadIndices.push_back(indexArray[cx+1][cy][cz]);
                if (cy-1 >= 0 && indexArray[cx][cy-1][cz] != -1) quadIndices.push_back(indexArray[cx][cy-1][cz]);
                if (cy+1 < hy && indexArray[cx][cy+1][cz] != -1) quadIndices.push_back(indexArray[cx][cy+1][cz]);
                if (cz-1 >= 0 && indexArray[cx][cy][cz-1] != -1) quadIndices.push_back(indexArray[cx][cy][cz-1]);
                if (cz+1 < hz && indexArray[cx][cy][cz+1] != -1) quadIndices.push_back(indexArray[cx][cy][cz+1]);
                
                result.cc[3].push_back(quadIndices);

                if (quadIndices.size() != 6) {
                    std::cout << "cube has " << quadIndices.size() << " quads (error)" << std::endl;
                    continue;
                }
                cubeIdx += 1;
                result.vals[3].push_back(val);
                
            }
        }
    }
    std::cout << "here line 779" << std::endl;
    // Clear indexArray to free memory (no longer needed after building cell complex)
    indexArray.clear();
    indexArray.shrink_to_fit();
    
    // Shrink result vectors to actual size
    for (int d = 0; d < 4; ++d) {
        result.cc[d].shrink_to_fit();
        result.vals[d].shrink_to_fit();
    }
    result.vertexCoords.shrink_to_fit();
    vertNeedMidPoint.shrink_to_fit();

    return {result, vertNeedMidPoint};
}

// -----------------------------------------------------------------------------
// Refinement of cubes into a mixed 3D complex (cubes + pyramids, optionally tets)
// -----------------------------------------------------------------------------

struct RefinementResult3D {
    std::vector<Point3D> vertexCoords;   // all vertices (original + new)
    std::vector<double> vertexVals;     // scalar values at refinedVerts
    std::vector<std::vector<int>> subs;           // top-dimensional cells (cubes/pyramids/tets)
};


// Refine cubes: simple cubes remain cubes; non-simple cubes are split into 6 pyramids.
// "Simplicity" criterion is intentionally simple here: val <= 0 means simple, otherwise refine.
RefinementResult3D refineCubesToMixed(const CellComplex& cc) {
    RefinementResult3D out;

    // Use references to avoid copying large structures (complex and vals are read-only)
    const std::vector<std::vector<std::vector<int>>>& complex = cc.cc;
    const std::vector<std::vector<double>>& vals = cc.vals;
    
    // These need to be copied since we modify them
    std::vector<Point3D> vertexCoords = cc.vertexCoords;
    std::vector<double> vertexVals = cc.vals[0];

    int tot = complex[0].size() + complex[1].size() + complex[2].size() + complex[3].size();
    vertexCoords.reserve(tot);

    int n = static_cast<int>(complex[0].size());
    int k = 4;  // dimensions: 0=vertices, 1=edges, 2=triangles, 3=tets

    std::vector<Point3D> ncoords(tot);
    std::vector<double> nvals(tot, 0.0);
    int ct = 0;  // new vertex counter

    auto getVinds = [&](int dim, int cellIdx, const std::vector<std::vector<std::vector<int>>>& currentSubs) -> std::vector<int> {
        if (dim == 1) {
            if (cellIdx < 0 || cellIdx >= static_cast<int>(complex[1].size())) {
                std::cerr << "Error in getVinds: cellIdx " << cellIdx << " out of bounds for cc[1].size()=" << complex[1].size() << std::endl;
                return {};
            }
            return cc.cc[1][cellIdx];  // edges are directly vertex indices
        } else {
            // For higher dimensions, get vertices from face subdivisions
            if (cellIdx < 0 || cellIdx >= static_cast<int>(complex[dim].size())) {
                std::cerr << "Error in getVinds: cellIdx " << cellIdx << " out of bounds for cc[" << dim << "].size()=" << complex[dim].size() << std::endl;
                return {};
            }
            const auto& faces = complex[dim][cellIdx];
            std::unordered_set<int> verts;
            for (int faceIdx : faces) {
                if (faceIdx < 0 || faceIdx >= static_cast<int>(currentSubs.size())) {
                    std::cerr << "Warning in getVinds: faceIdx " << faceIdx << " out of bounds for currentSubs.size()=" << currentSubs.size() << ", dim=" << dim << ", cellIdx=" << cellIdx << std::endl;
                    continue;
                }
                for (const auto& sub : currentSubs[faceIdx]) {
                    for (int v : sub) {
                        verts.insert(v);
                    }
                }
            }
            return std::vector<int>(verts.begin(), verts.end());
        }
    };


    // Initialize subdivisions: for each original cell, store list of refined cells (as vertex indices)
    // Start with vertex subdivisions (each vertex is itself)
    std::vector<std::vector<std::vector<int>>> subs(vertexCoords.size());
    for (size_t i = 0; i < vertexCoords.size(); ++i) {
        subs[i] = {{static_cast<int>(i)}};
    }


    for (int i = 1; i < k; ++i) {
        int dim = i;
        std::vector<std::vector<std::vector<int>>> nsubs(complex[dim].size());
        std::cout << "Processing dim=" << dim << ", cc[" << dim << "].size()=" << complex[dim].size() 
                  << ", subs.size()=" << subs.size() << std::endl;

        for (size_t j = 0; j < complex[dim].size(); ++j) {

            const auto& faces = complex[dim][j];  // face indices (lower dimension)
            if (j >= vals[dim].size()) {
                std::cerr << "Error: j=" << j << " out of bounds for vals[" << dim << "].size()=" << vals[dim].size() << std::endl;
                continue;
            }
            double val = vals[dim][j];  // cell value

            // Get face values
            std::vector<double> facevals;
            facevals.reserve(faces.size());
            for (int faceIdx : faces) {
                if (faceIdx < 0 || faceIdx >= static_cast<int>(vals[dim-1].size())) {
                    std::cerr << "Error: faceIdx " << faceIdx << " out of bounds for vals[" << (dim-1) 
                              << "].size()=" << vals[dim-1].size() << ", dim=" << dim << ", j=" << j << std::endl;
                    continue;
                }
                facevals.push_back(vals[dim-1][faceIdx]);
            }
            facevals.shrink_to_fit();
            
            // Get face subdivisions (from current subs, which contains subdivisions of dimension dim-1)
            std::vector<std::vector<std::vector<int>>> facesubs;
            facesubs.reserve(faces.size()*3);
            for (int faceIdx : faces) {
                if (faceIdx < 0 || faceIdx >= static_cast<int>(subs.size())) {
                    std::cerr << "Error: faceIdx " << faceIdx << " out of bounds for dim=" << dim 
                              << ", subs.size()=" << subs.size() << ", j=" << j << std::endl;
                    continue;
                }
                facesubs.push_back(subs[faceIdx]);
            }
            facesubs.shrink_to_fit();
            // Get vertex indices for this cell using getVinds helper
            std::vector<int> verts = getVinds(dim, j, subs);

            if (verts.empty()) {
                std::cerr << "Warning: empty verts for dim=" << dim << ", j=" << j << std::endl;
                nsubs[j] = {};
                continue;
            }


            bool isSimple = true;
            if (val <= 0.0) {
                isSimple = true;
            } else {
                if (!facevals.empty()) {
                    // Check if Max[facevals] > 0
                    double maxFaceVal = *std::max_element(facevals.begin(), facevals.end());
                    if (maxFaceVal > 0.0) {

                        int maxSubLength = 0;
                        for (const auto& faceSub : facesubs) {
                            maxSubLength = std::max(maxSubLength, static_cast<int>(faceSub.size()));
                        }
                        if (maxSubLength == 1) {
                            isSimple = true;
                        }
                        else{
                            isSimple = false;
                        }
                        
                    }
                    else{
                        isSimple = false;
                    }
                    
                }
                else{
                    std::cout <<"Facevals empty for cell (this shouldnt print)"<< j << std::endl;
                }
            }
            
            if (isSimple) {
                // Cell is simple: use original cell (verts)
                nsubs[j] = {verts};
            }
            else {
                // Cell is not simple: create new vertex, and form one cell connecting the vertex to each face subdivision
                // Compute centroid
                Point3D centroid(0.0, 0.0, 0.0);
                
                for (int vIdx : verts) {
                    if (vIdx < static_cast<int>(vertexCoords.size())) {
                        // Original vertex
                        centroid.x += vertexCoords[vIdx].x;
                        centroid.y += vertexCoords[vIdx].y;
                        centroid.z += vertexCoords[vIdx].z;
                    } else {
                        // New refined vertex (from previous refinement steps)
                        int ncoordsIdx = vIdx - static_cast<int>(vertexCoords.size());
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
                centroid.x /= verts.size();
                centroid.y /= verts.size();
                centroid.z /= verts.size();

                if (ct >= static_cast<int>(ncoords.size())) {
                    std::cerr << "Error: ct=" << ct << " >= ncoords.size()=" << ncoords.size() 
                              << ", dim=" << dim << ", j=" << j << std::endl;
                    throw std::runtime_error("New vertex counter out of bounds");
                }
                
                ncoords[ct] = centroid;
                nvals[ct] = val;
                int newVertIdx = n + ct;
                ct += 1;
                
                // Form new cells: connect new vertex to each face subdivision
                std::vector<std::vector<int>> newCells;
                newCells.reserve(facesubs.size()*5);
                for (const auto& faceSub : facesubs) {
                    for (const auto& sub : faceSub) {
                        if (sub.empty()) {
                            std::cerr << "Warning: empty subdivision for dim=" << dim << ", j=" << j << std::endl;
                            continue;
                        }
                        std::vector<int> newCell = sub;
                        
                        newCell.push_back(newVertIdx);
                        // For tets (dim=3), each new cell should have 4 vertices
                        if (dim == 3 && (newCell.size() != 4 && newCell.size() != 5 && newCell.size() != 8)) {
                            std::cerr << "Warning: tet cell has " << newCell.size() 
                                      << " vertices instead of 4/5/8 for dim=" << dim << ", j=" << j << std::endl;
                        }

                        newCells.push_back(newCell);
                    }
                }
                newCells.shrink_to_fit();

                nsubs[j] = newCells;
                
            }
        }
        
        // Update subs for next iteration (Mathematica: subs = nsubs)
        subs = std::move(nsubs);  // Use move to avoid copying
        subs.shrink_to_fit();
        std::cout<<"ct: "<< ct << std::endl;
        

    }

    //subs is all tets

    vertexCoords.insert(vertexCoords.end(), ncoords.begin(), ncoords.begin() + ct);
    vertexVals.insert(vertexVals.end(), nvals.begin(), nvals.begin() + ct);
    
    // Free intermediate arrays
    ncoords.clear();
    nvals.clear();
    ncoords.shrink_to_fit();
    nvals.shrink_to_fit();

    std::vector<std::vector<int>> threeDimensionalCells;
    threeDimensionalCells.reserve(subs.size()*2);

    for (const auto& tetSub : subs) {  // subs now contains tet subdivisions
        for (const auto& sub : tetSub) {

            if (sub.size() == 4){
                threeDimensionalCells.push_back(sub);
            }
            else if (sub.size() == 8){

                std::vector<int> cubeVerts = sub;

                // Compute cube center in high-res coordinates
                double minX = std::numeric_limits<int>::max();
                double minY = std::numeric_limits<int>::max();
                double minZ = std::numeric_limits<int>::max();
                double maxX = std::numeric_limits<int>::lowest();
                double maxY = std::numeric_limits<int>::lowest();
                double maxZ = std::numeric_limits<int>::lowest();

                for (int v : sub) {
                    double vx = vertexCoords[v].x;
                    double vy = vertexCoords[v].y;
                    double vz = vertexCoords[v].z;

                    minX = std::min(minX, vx);
                    minY = std::min(minY, vy);
                    minZ = std::min(minZ, vz);
                    maxX = std::max(maxX, vx);
                    maxY = std::max(maxY, vy);
                    maxZ = std::max(maxZ, vz);
                }

                double cx = (minX + maxX) / 2;
                double cy = (minY + maxY) / 2;
                double cz = (minZ + maxZ) / 2;

                // Assign each vertex to a slot based on its position relative to center.
                // The slot index must match the offset order in marching_cubes.cpp:
                // (dx,dy,dz) in {0,1}^3 is mapped as:
                //  (0,0,0)->0, (0,0,1)->1, (0,1,0)->2, (0,1,1)->3,
                //  (1,0,0)->4, (1,0,1)->5, (1,1,0)->6, (1,1,1)->7.
                for (int v : sub) {
                    double vx = vertexCoords[v].x;
                    double vy = vertexCoords[v].y;
                    double vz = vertexCoords[v].z;
                    int dx = (vx > cx) ? 1 : 0;
                    int dy = (vy > cy) ? 1 : 0;
                    int dz = (vz > cz) ? 1 : 0;
                    int slot = dx * 4 + dy * 2 + dz; // 0..7, consistent with offset table
                    if (slot < 0 || slot >= 8) continue;
                    cubeVerts[slot] = v;
                }

                threeDimensionalCells.push_back(cubeVerts);
            }
            else if (sub.size() == 5){
                // Canonicalize pyramid vertex ordering:
                // local index 0 = apex, 1..4 = base quad in pattern:
                //  (0,0,0)->1, (0,0,1)->2, (0,1,0)->3, (0,1,1)->4
                std::vector<int> pyr(5, -1);
                const double eps = 1e-6;

                // Compute bounding box
                double minX = std::numeric_limits<double>::max();
                double minY = std::numeric_limits<double>::max();
                double minZ = std::numeric_limits<double>::max();
                double maxX = -std::numeric_limits<double>::max();
                double maxY = -std::numeric_limits<double>::max();
                double maxZ = -std::numeric_limits<double>::max();

                for (int v : sub) {
                    double vx = vertexCoords[v].x;
                    double vy = vertexCoords[v].y;
                    double vz = vertexCoords[v].z;
                    minX = std::min(minX, vx);
                    minY = std::min(minY, vy);
                    minZ = std::min(minZ, vz);
                    maxX = std::max(maxX, vx);
                    maxY = std::max(maxY, vy);
                    maxZ = std::max(maxZ, vz);
                }

                // Identify apex: the vertex that is not on the bounding-box corners
                int apex = -1;
                std::vector<int> base;
                base.reserve(4);
                for (int v : sub) {
                    double vx = vertexCoords[v].x;
                    double vy = vertexCoords[v].y;
                    double vz = vertexCoords[v].z;
                    bool onX = (std::abs(vx - minX) < eps) || (std::abs(vx - maxX) < eps);
                    bool onY = (std::abs(vy - minY) < eps) || (std::abs(vy - maxY) < eps);
                    bool onZ = (std::abs(vz - minZ) < eps) || (std::abs(vz - maxZ) < eps);
                    if (onX && onY && onZ) {
                        base.push_back(v);
                    } else {
                        apex = v;
                    }
                }
                if (apex == -1 || base.size() != 4) {
                    // Fallback: keep original order
                    threeDimensionalCells.push_back(sub);
                    std::cout << "NO APEX!! This should not print, apex is " << apex << "base size is " << base.size() << std::endl;

                    continue;
                }

                pyr[0] = apex;

                // Determine which coordinate is constant across the base (the pyramid's normal)
                double bx0 = vertexCoords[base[0]].x;
                double by0 = vertexCoords[base[0]].y;
                double bz0 = vertexCoords[base[0]].z;

                int constDim = -1;
                if (std::all_of(base.begin(), base.end(), [&](int v){
                        return std::abs(vertexCoords[v].x - bx0) < eps;
                    })) {
                    constDim = 0; // x constant
                } else if (std::all_of(base.begin(), base.end(), [&](int v){
                        return std::abs(vertexCoords[v].y - by0) < eps;
                    })) {
                    constDim = 1; // y constant
                } else if (std::all_of(base.begin(), base.end(), [&](int v){
                        return std::abs(vertexCoords[v].z - bz0) < eps;
                    })) {
                    constDim = 2; // z constant
                } else {
                    // Degenerate; keep as-is
                    threeDimensionalCells.push_back(sub);
                    continue;
                }

                // Choose the two varying dims as (u,v)
                int uDim = (constDim == 0) ? 1 : 0;
                int vDim = (constDim == 2) ? 1 : 2;
                if (constDim == 1) {
                    uDim = 0;
                    vDim = 2;
                }

                // Compute min/max for u,v over base
                double uMin = std::numeric_limits<double>::max();
                double vMin = std::numeric_limits<double>::max();
                double uMax = -std::numeric_limits<double>::max();
                double vMax = -std::numeric_limits<double>::max();

                auto getCoord = [&](int v, int dim){
                    if (dim == 0) return vertexCoords[v].x;
                    if (dim == 1) return vertexCoords[v].y;
                    return vertexCoords[v].z;
                };

                for (int v : base) {
                    double u = getCoord(v, uDim);
                    double vv = getCoord(v, vDim);
                    uMin = std::min(uMin, u);
                    vMin = std::min(vMin, vv);
                    uMax = std::max(uMax, u);
                    vMax = std::max(vMax, vv);
                }

                // Assign the four base vertices to slots 1..4 according to (u,v):
                // (0,0)->1, (0,1)->2, (1,0)->3, (1,1)->4
                for (int v : base) {
                    double u = getCoord(v, uDim);
                    double vv = getCoord(v, vDim);
                    int uBit = (std::abs(u - uMin) < eps) ? 0 : 1;
                    int vBit = (std::abs(vv - vMin) < eps) ? 0 : 1;
                    int slot = -1;
                    if (uBit == 0 && vBit == 0) slot = 1;
                    else if (uBit == 0 && vBit == 1) slot = 2;
                    else if (uBit == 1 && vBit == 0) slot = 3;
                    else if (uBit == 1 && vBit == 1) slot = 4;
                    if (slot >= 1 && slot <= 4) {
                        pyr[slot] = v;
                    }
                }

                // If any slot not assigned, fallback to original order
                bool ok = true;
                for (int i = 0; i < 5; ++i) {
                    if (pyr[i] == -1) { ok = false; break; }
                }
                if (!ok) {
                    std::cout << "pyramid processing not OK!! This should not print";
                    threeDimensionalCells.push_back(sub);
                } else {
                    threeDimensionalCells.push_back(pyr);
                }
            }

        }
    }

    // Clear subs to free memory before returning
    subs.clear();
    subs.shrink_to_fit();
    
    out.vertexCoords = std::move(vertexCoords);
    out.vertexVals = std::move(vertexVals);
    out.subs = std::move(threeDimensionalCells);
    
    // Shrink to fit final result
    out.vertexCoords.shrink_to_fit();
    out.vertexVals.shrink_to_fit();
    out.subs.shrink_to_fit();
    
    return out;
}

// -----------------------------------------------------------------------------
// Mixed marching routine on a refined 3D cell complex
// -----------------------------------------------------------------------------
//
// We assume three separate lookup tables (to be defined elsewhere):
//  - cubeCycles   : for 8-vertex cubes
//  - tetCycles    : for 4-vertex tetrahedra
//  - pyramidCycles: for 5-vertex pyramids
//
// Each table is indexed by a bit key built from the sign of the scalar values
// at the ordered vertices of the corresponding cell type.
//
// The local vertex ordering for each cell in Cell3D::verts must match the
// convention assumed by the corresponding lookup table.
// -----------------------------------------------------------------------------

// These will be defined and filled elsewhere (e.g. in marching_cubes.cpp or a new file).


std::pair<std::vector<Point3D>, std::vector<Triangle>>
contourMixedCells(const std::vector<Point3D>& vts,
                  const std::vector<double>& vals,
                  std::vector<std::vector<int>> cells,
                  int originalNumVertices,
                  std::vector<bool>vertNeedMidPoint
                )
{   

    std::vector<Point3D>& zeroVerts = std::vector<Point3D>();
    zeroVerts.reserve(vts.size());

    std::vector<std::pair<Point3D, Point3D>> zeroCrossingEdges;
    zeroCrossingEdges.reserve(vts.size());

    std::vector<bool>usedMidPoint;
    usedMidPoint.reserve(vts.size());

    std::vector<Triangle>& zeroTris = std::vector<Triangle>();
    zeroTris.reserve(cells.size());
    
    std::unordered_map<UndirectedEdgeKey, int, UndirectedEdgeKeyHash> edgeHash;

    std::vector<std::vector<std::vector<std::pair<int,int>>>> cubeCycles(256);

    cubeCycles = {{}, {{{4, 8}, {6, 8}, {7, 8}}}, {{{3, 7}, {7, 8}, {5, 7}}}, 
    {{{5, 7}, {3, 7}, {4, 8}}, {{6, 8}, {5, 7}, {4, 8}}}, 
    {{{2, 6}, {5, 6}, {6, 8}}}, {{{2, 6}, {5, 6}, {7, 8}}, 
     {{4, 8}, {2, 6}, {7, 8}}}, {{{2, 6}, {5, 6}, {5, 7}}, 
     {{3, 7}, {2, 6}, {5, 7}}, {{6, 8}, {2, 6}, {3, 7}}, 
     {{7, 8}, {6, 8}, {3, 7}}}, {{{2, 6}, {3, 7}, {4, 8}}, 
     {{2, 6}, {5, 6}, {5, 7}}, {{3, 7}, {2, 6}, {5, 7}}}, 
    {{{1, 5}, {5, 7}, {5, 6}}}, {{{1, 5}, {5, 7}, {7, 8}}, 
     {{4, 8}, {1, 5}, {7, 8}}, {{5, 6}, {1, 5}, {4, 8}}, 
     {{6, 8}, {5, 6}, {4, 8}}}, {{{3, 7}, {7, 8}, {5, 6}}, 
     {{1, 5}, {3, 7}, {5, 6}}}, {{{1, 5}, {3, 7}, {4, 8}}, 
     {{4, 8}, {6, 8}, {5, 6}}, {{1, 5}, {4, 8}, {5, 6}}}, 
    {{{1, 5}, {5, 7}, {6, 8}}, {{2, 6}, {1, 5}, {6, 8}}}, 
    {{{2, 6}, {1, 5}, {4, 8}}, {{1, 5}, {5, 7}, {7, 8}}, 
     {{4, 8}, {1, 5}, {7, 8}}}, {{{2, 6}, {1, 5}, {3, 7}}, 
     {{3, 7}, {7, 8}, {6, 8}}, {{2, 6}, {3, 7}, {6, 8}}}, 
    {{{4, 8}, {2, 6}, {1, 5}}, {{3, 7}, {4, 8}, {1, 5}}}, 
    {{{4, 8}, {3, 4}, {2, 4}}}, {{{3, 4}, {2, 4}, {6, 8}}, 
     {{7, 8}, {3, 4}, {6, 8}}}, {{{2, 4}, {4, 8}, {7, 8}}, 
     {{5, 7}, {2, 4}, {7, 8}}, {{3, 4}, {2, 4}, {5, 7}}, 
     {{3, 7}, {3, 4}, {5, 7}}}, {{{2, 4}, {6, 8}, {5, 7}}, 
     {{5, 7}, {3, 7}, {3, 4}}, {{2, 4}, {5, 7}, {3, 4}}}, 
    {{{5, 6}, {6, 8}, {4, 8}}, {{3, 4}, {5, 6}, {4, 8}}, 
     {{2, 6}, {5, 6}, {3, 4}}, {{2, 4}, {2, 6}, {3, 4}}}, 
    {{{3, 4}, {5, 6}, {7, 8}}, {{3, 4}, {2, 4}, {2, 6}}, 
     {{5, 6}, {3, 4}, {2, 6}}}, {{{6, 8}, {4, 8}, {7, 8}}, 
     {{3, 4}, {2, 4}, {2, 6}}, {{3, 7}, {3, 4}, {2, 6}}, 
     {{5, 7}, {3, 7}, {2, 6}}, {{5, 6}, {5, 7}, {2, 6}}}, 
    {{{3, 4}, {2, 4}, {2, 6}}, {{3, 7}, {3, 4}, {2, 6}}, 
     {{5, 7}, {3, 7}, {2, 6}}, {{5, 6}, {5, 7}, {2, 6}}}, 
    {{{2, 4}, {1, 5}, {3, 4}}, {{2, 4}, {5, 6}, {1, 5}}, 
     {{2, 4}, {4, 8}, {5, 6}}, {{5, 6}, {4, 8}, {5, 7}}, 
     {{3, 4}, {5, 7}, {4, 8}}, {{3, 4}, {1, 5}, {5, 7}}}, 
    {{{2, 4}, {1, 5}, {3, 4}}, {{2, 4}, {5, 6}, {1, 5}}, 
     {{2, 4}, {6, 8}, {5, 6}}, {{3, 4}, {5, 7}, {7, 8}}, 
     {{3, 4}, {1, 5}, {5, 7}}}, {{{3, 4}, {2, 4}, {1, 5}}, 
     {{3, 4}, {1, 5}, {3, 7}}, {{5, 6}, {1, 5}, {2, 4}}, 
     {{5, 6}, {4, 8}, {7, 8}}, {{5, 6}, {2, 4}, {4, 8}}}, 
    {{{3, 4}, {2, 4}, {1, 5}}, {{2, 4}, {5, 6}, {1, 5}}, 
     {{3, 4}, {1, 5}, {3, 7}}, {{2, 4}, {6, 8}, {5, 6}}}, 
    {{{2, 4}, {1, 5}, {3, 4}}, {{2, 4}, {2, 6}, {1, 5}}, 
     {{3, 4}, {1, 5}, {5, 7}}, {{6, 8}, {4, 8}, {5, 7}}, 
     {{3, 4}, {5, 7}, {4, 8}}}, {{{3, 4}, {5, 7}, {7, 8}}, 
     {{2, 4}, {2, 6}, {1, 5}}, {{3, 4}, {2, 4}, {1, 5}}, 
     {{3, 4}, {1, 5}, {5, 7}}}, {{{3, 4}, {1, 5}, {3, 7}}, 
     {{6, 8}, {4, 8}, {7, 8}}, {{3, 4}, {2, 4}, {1, 5}}, 
     {{1, 5}, {2, 4}, {2, 6}}}, {{{3, 4}, {1, 5}, {3, 7}}, 
     {{3, 4}, {2, 4}, {1, 5}}, {{1, 5}, {2, 4}, {2, 6}}}, 
    {{{3, 7}, {1, 3}, {3, 4}}}, {{{6, 8}, {7, 8}, {3, 7}}, 
     {{1, 3}, {6, 8}, {3, 7}}, {{4, 8}, {6, 8}, {1, 3}}, 
     {{3, 4}, {4, 8}, {1, 3}}}, {{{1, 3}, {3, 4}, {7, 8}}, 
     {{5, 7}, {1, 3}, {7, 8}}}, {{{1, 3}, {6, 8}, {5, 7}}, 
     {{1, 3}, {3, 4}, {4, 8}}, {{6, 8}, {1, 3}, {4, 8}}}, 
    {{{1, 3}, {3, 4}, {2, 6}}, {{1, 3}, {2, 6}, {5, 6}}, 
     {{1, 3}, {5, 6}, {3, 7}}, {{5, 6}, {6, 8}, {3, 7}}, 
     {{3, 4}, {3, 7}, {6, 8}}, {{3, 4}, {6, 8}, {2, 6}}}, 
    {{{3, 4}, {2, 6}, {1, 3}}, {{3, 4}, {4, 8}, {2, 6}}, 
     {{5, 6}, {1, 3}, {2, 6}}, {{5, 6}, {7, 8}, {3, 7}}, 
     {{5, 6}, {3, 7}, {1, 3}}}, {{{3, 4}, {2, 6}, {1, 3}}, 
     {{3, 4}, {6, 8}, {2, 6}}, {{1, 3}, {2, 6}, {5, 6}}, 
     {{3, 4}, {7, 8}, {6, 8}}, {{1, 3}, {5, 6}, {5, 7}}}, 
    {{{1, 3}, {5, 6}, {5, 7}}, {{1, 3}, {3, 4}, {2, 6}}, 
     {{1, 3}, {2, 6}, {5, 6}}, {{3, 4}, {4, 8}, {2, 6}}}, 
    {{{3, 4}, {3, 7}, {5, 7}}, {{5, 6}, {3, 4}, {5, 7}}, 
     {{1, 3}, {3, 4}, {5, 6}}, {{1, 5}, {1, 3}, {5, 6}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}, {{4, 8}, {6, 8}, {5, 6}}, 
     {{3, 4}, {4, 8}, {5, 6}}, {{1, 3}, {3, 4}, {5, 6}}, 
     {{1, 5}, {1, 3}, {5, 6}}}, {{{3, 4}, {7, 8}, {5, 6}}, 
     {{5, 6}, {1, 5}, {1, 3}}, {{3, 4}, {5, 6}, {1, 3}}}, 
    {{{4, 8}, {6, 8}, {5, 6}}, {{3, 4}, {4, 8}, {5, 6}}, 
     {{1, 3}, {3, 4}, {5, 6}}, {{1, 5}, {1, 3}, {5, 6}}}, 
    {{{1, 3}, {3, 4}, {2, 6}}, {{1, 3}, {2, 6}, {1, 5}}, 
     {{6, 8}, {2, 6}, {3, 4}}, {{6, 8}, {3, 7}, {5, 7}}, 
     {{6, 8}, {3, 4}, {3, 7}}}, {{{1, 3}, {2, 6}, {1, 5}}, 
     {{7, 8}, {3, 7}, {5, 7}}, {{1, 3}, {3, 4}, {2, 6}}, 
     {{2, 6}, {3, 4}, {4, 8}}}, {{{1, 3}, {3, 4}, {2, 6}}, 
     {{3, 4}, {6, 8}, {2, 6}}, {{1, 3}, {2, 6}, {1, 5}}, 
     {{3, 4}, {7, 8}, {6, 8}}}, {{{1, 3}, {2, 6}, {1, 5}}, 
     {{1, 3}, {3, 4}, {2, 6}}, {{2, 6}, {3, 4}, {4, 8}}}, 
    {{{2, 4}, {4, 8}, {3, 7}}, {{1, 3}, {2, 4}, {3, 7}}}, 
    {{{2, 4}, {6, 8}, {1, 3}}, {{6, 8}, {7, 8}, {3, 7}}, 
     {{1, 3}, {6, 8}, {3, 7}}}, {{{2, 4}, {5, 7}, {1, 3}}, 
     {{2, 4}, {4, 8}, {7, 8}}, {{5, 7}, {2, 4}, {7, 8}}}, 
    {{{5, 7}, {1, 3}, {2, 4}}, {{6, 8}, {5, 7}, {2, 4}}}, 
    {{{2, 4}, {2, 6}, {1, 3}}, {{1, 3}, {2, 6}, {5, 6}}, 
     {{6, 8}, {4, 8}, {3, 7}}, {{6, 8}, {3, 7}, {5, 6}}, 
     {{1, 3}, {5, 6}, {3, 7}}}, {{{5, 6}, {7, 8}, {3, 7}}, 
     {{5, 6}, {3, 7}, {1, 3}}, {{2, 4}, {2, 6}, {1, 3}}, 
     {{5, 6}, {1, 3}, {2, 6}}}, {{{1, 3}, {2, 6}, {5, 6}}, 
     {{1, 3}, {5, 6}, {5, 7}}, {{1, 3}, {2, 4}, {2, 6}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{1, 3}, {2, 6}, {5, 6}}, 
     {{1, 3}, {5, 6}, {5, 7}}, {{1, 3}, {2, 4}, {2, 6}}}, 
    {{{2, 4}, {4, 8}, {5, 6}}, {{5, 7}, {5, 6}, {4, 8}}, 
     {{5, 7}, {4, 8}, {3, 7}}, {{2, 4}, {1, 5}, {1, 3}}, 
     {{2, 4}, {5, 6}, {1, 5}}}, {{{2, 4}, {1, 5}, {1, 3}}, 
     {{7, 8}, {3, 7}, {5, 7}}, {{2, 4}, {6, 8}, {5, 6}}, 
     {{2, 4}, {5, 6}, {1, 5}}}, {{{2, 4}, {1, 5}, {1, 3}}, 
     {{2, 4}, {5, 6}, {1, 5}}, {{2, 4}, {4, 8}, {5, 6}}, 
     {{5, 6}, {4, 8}, {7, 8}}}, {{{2, 4}, {1, 5}, {1, 3}}, 
     {{2, 4}, {6, 8}, {5, 6}}, {{2, 4}, {5, 6}, {1, 5}}}, 
    {{{1, 5}, {1, 3}, {2, 4}}, {{2, 6}, {1, 5}, {2, 4}}, 
     {{5, 7}, {6, 8}, {4, 8}}, {{3, 7}, {5, 7}, {4, 8}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}, {{1, 5}, {1, 3}, {2, 4}}, 
     {{2, 6}, {1, 5}, {2, 4}}}, {{{6, 8}, {4, 8}, {7, 8}}, 
     {{1, 5}, {1, 3}, {2, 4}}, {{2, 6}, {1, 5}, {2, 4}}}, 
    {{{1, 5}, {1, 3}, {2, 4}}, {{2, 6}, {1, 5}, {2, 4}}}, 
    {{{2, 6}, {2, 4}, {1, 2}}}, {{{1, 2}, {2, 6}, {6, 8}}, 
     {{7, 8}, {1, 2}, {6, 8}}, {{2, 4}, {1, 2}, {7, 8}}, 
     {{4, 8}, {2, 4}, {7, 8}}}, {{{2, 4}, {1, 2}, {3, 7}}, 
     {{2, 4}, {3, 7}, {7, 8}}, {{2, 4}, {7, 8}, {2, 6}}, 
     {{7, 8}, {5, 7}, {2, 6}}, {{1, 2}, {2, 6}, {5, 7}}, 
     {{1, 2}, {5, 7}, {3, 7}}}, {{{2, 4}, {1, 2}, {3, 7}}, 
     {{2, 4}, {3, 7}, {4, 8}}, {{1, 2}, {5, 7}, {3, 7}}, 
     {{6, 8}, {5, 7}, {2, 6}}, {{1, 2}, {2, 6}, {5, 7}}}, 
    {{{5, 6}, {6, 8}, {2, 4}}, {{1, 2}, {5, 6}, {2, 4}}}, 
    {{{1, 2}, {5, 6}, {7, 8}}, {{7, 8}, {4, 8}, {2, 4}}, 
     {{1, 2}, {7, 8}, {2, 4}}}, {{{1, 2}, {3, 7}, {2, 4}}, 
     {{1, 2}, {5, 6}, {5, 7}}, {{1, 2}, {5, 7}, {3, 7}}, 
     {{2, 4}, {7, 8}, {6, 8}}, {{2, 4}, {3, 7}, {7, 8}}}, 
    {{{1, 2}, {5, 6}, {5, 7}}, {{2, 4}, {3, 7}, {4, 8}}, 
     {{1, 2}, {3, 7}, {2, 4}}, {{1, 2}, {5, 7}, {3, 7}}}, 
    {{{5, 7}, {5, 6}, {2, 6}}, {{2, 4}, {5, 7}, {2, 6}}, 
     {{1, 5}, {5, 7}, {2, 4}}, {{1, 2}, {1, 5}, {2, 4}}}, 
    {{{5, 6}, {2, 6}, {6, 8}}, {{1, 5}, {5, 7}, {7, 8}}, 
     {{1, 2}, {1, 5}, {7, 8}}, {{2, 4}, {1, 2}, {7, 8}}, 
     {{4, 8}, {2, 4}, {7, 8}}}, {{{1, 2}, {3, 7}, {2, 4}}, 
     {{1, 2}, {1, 5}, {3, 7}}, {{5, 6}, {2, 6}, {7, 8}}, 
     {{2, 4}, {7, 8}, {2, 6}}, {{2, 4}, {3, 7}, {7, 8}}}, 
    {{{1, 2}, {1, 5}, {3, 7}}, {{1, 2}, {3, 7}, {2, 4}}, 
     {{5, 6}, {2, 6}, {6, 8}}, {{2, 4}, {3, 7}, {4, 8}}}, 
    {{{2, 4}, {5, 7}, {6, 8}}, {{2, 4}, {1, 2}, {1, 5}}, 
     {{5, 7}, {2, 4}, {1, 5}}}, {{{7, 8}, {4, 8}, {2, 4}}, 
     {{5, 7}, {7, 8}, {2, 4}}, {{1, 5}, {5, 7}, {2, 4}}, 
     {{1, 2}, {1, 5}, {2, 4}}}, {{{1, 2}, {1, 5}, {3, 7}}, 
     {{1, 2}, {3, 7}, {2, 4}}, {{2, 4}, {7, 8}, {6, 8}}, 
     {{2, 4}, {3, 7}, {7, 8}}}, {{{1, 2}, {1, 5}, {3, 7}}, 
     {{1, 2}, {3, 7}, {2, 4}}, {{2, 4}, {3, 7}, {4, 8}}}, 
    {{{1, 2}, {2, 6}, {4, 8}}, {{3, 4}, {1, 2}, {4, 8}}}, 
    {{{1, 2}, {7, 8}, {3, 4}}, {{1, 2}, {2, 6}, {6, 8}}, 
     {{7, 8}, {1, 2}, {6, 8}}}, {{{7, 8}, {5, 7}, {2, 6}}, 
     {{7, 8}, {2, 6}, {4, 8}}, {{1, 2}, {2, 6}, {5, 7}}, 
     {{1, 2}, {3, 7}, {3, 4}}, {{1, 2}, {5, 7}, {3, 7}}}, 
    {{{1, 2}, {3, 7}, {3, 4}}, {{1, 2}, {5, 7}, {3, 7}}, 
     {{6, 8}, {5, 7}, {2, 6}}, {{1, 2}, {2, 6}, {5, 7}}}, 
    {{{1, 2}, {5, 6}, {3, 4}}, {{5, 6}, {6, 8}, {4, 8}}, 
     {{3, 4}, {5, 6}, {4, 8}}}, {{{3, 4}, {1, 2}, {5, 6}}, 
     {{7, 8}, {3, 4}, {5, 6}}}, {{{1, 2}, {3, 7}, {3, 4}}, 
     {{1, 2}, {5, 6}, {5, 7}}, {{1, 2}, {5, 7}, {3, 7}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{1, 2}, {3, 7}, {3, 4}}, 
     {{1, 2}, {5, 6}, {5, 7}}, {{1, 2}, {5, 7}, {3, 7}}}, 
    {{{1, 2}, {1, 5}, {3, 4}}, {{3, 4}, {1, 5}, {5, 7}}, 
     {{5, 6}, {2, 6}, {4, 8}}, {{5, 6}, {4, 8}, {5, 7}}, 
     {{3, 4}, {5, 7}, {4, 8}}}, {{{5, 6}, {2, 6}, {6, 8}}, 
     {{1, 2}, {1, 5}, {3, 4}}, {{7, 8}, {3, 4}, {5, 7}}, 
     {{3, 4}, {1, 5}, {5, 7}}}, {{{3, 4}, {1, 2}, {1, 5}}, 
     {{3, 7}, {3, 4}, {1, 5}}, {{7, 8}, {5, 6}, {2, 6}}, 
     {{4, 8}, {7, 8}, {2, 6}}}, {{{5, 6}, {2, 6}, {6, 8}}, 
     {{3, 7}, {3, 4}, {1, 2}}, {{1, 5}, {3, 7}, {1, 2}}}, 
    {{{1, 2}, {1, 5}, {3, 4}}, {{6, 8}, {4, 8}, {5, 7}}, 
     {{3, 4}, {5, 7}, {4, 8}}, {{3, 4}, {1, 5}, {5, 7}}}, 
    {{{1, 2}, {1, 5}, {3, 4}}, {{7, 8}, {3, 4}, {5, 7}}, 
     {{3, 4}, {1, 5}, {5, 7}}}, {{{6, 8}, {4, 8}, {7, 8}}, 
     {{3, 7}, {3, 4}, {1, 2}}, {{1, 5}, {3, 7}, {1, 2}}}, 
    {{{3, 7}, {3, 4}, {1, 2}}, {{1, 5}, {3, 7}, {1, 2}}}, 
    {{{3, 7}, {1, 3}, {1, 2}}, {{2, 6}, {3, 7}, {1, 2}}, 
     {{3, 4}, {3, 7}, {2, 6}}, {{2, 4}, {3, 4}, {2, 6}}}, 
    {{{2, 4}, {3, 4}, {4, 8}}, {{1, 2}, {2, 6}, {6, 8}}, 
     {{1, 3}, {1, 2}, {6, 8}}, {{3, 7}, {1, 3}, {6, 8}}, 
     {{7, 8}, {3, 7}, {6, 8}}}, {{{1, 2}, {2, 6}, {5, 7}}, 
     {{1, 2}, {5, 7}, {1, 3}}, {{2, 4}, {3, 4}, {7, 8}}, 
     {{2, 4}, {7, 8}, {2, 6}}, {{5, 7}, {2, 6}, {7, 8}}}, 
    {{{1, 2}, {5, 7}, {1, 3}}, {{2, 4}, {3, 4}, {4, 8}}, 
     {{6, 8}, {5, 7}, {2, 6}}, {{1, 2}, {2, 6}, {5, 7}}}, 
    {{{1, 2}, {5, 6}, {1, 3}}, {{5, 6}, {3, 7}, {1, 3}}, 
     {{2, 4}, {3, 4}, {6, 8}}, {{3, 7}, {6, 8}, {3, 4}}, 
     {{5, 6}, {6, 8}, {3, 7}}}, {{{2, 4}, {3, 4}, {4, 8}}, 
     {{1, 2}, {5, 6}, {1, 3}}, {{7, 8}, {3, 7}, {5, 6}}, 
     {{5, 6}, {3, 7}, {1, 3}}}, {{{1, 3}, {1, 2}, {5, 6}}, 
     {{5, 7}, {1, 3}, {5, 6}}, {{6, 8}, {2, 4}, {3, 4}}, 
     {{7, 8}, {6, 8}, {3, 4}}}, {{{2, 4}, {3, 4}, {4, 8}}, 
     {{5, 7}, {1, 3}, {1, 2}}, {{5, 6}, {5, 7}, {1, 2}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{3, 4}, {3, 7}, {5, 7}}, 
     {{2, 4}, {3, 4}, {5, 7}}, {{2, 6}, {2, 4}, {5, 7}}, 
     {{5, 6}, {2, 6}, {5, 7}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{5, 6}, {2, 6}, {6, 8}}, {{2, 4}, {3, 4}, {4, 8}}, 
     {{7, 8}, {3, 7}, {5, 7}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{5, 6}, {2, 6}, {7, 8}}, {{2, 4}, {7, 8}, {2, 6}}, 
     {{2, 4}, {3, 4}, {7, 8}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{5, 6}, {2, 6}, {6, 8}}, {{3, 4}, {4, 8}, {2, 4}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{2, 4}, {3, 4}, {6, 8}}, 
     {{5, 7}, {6, 8}, {3, 7}}, {{6, 8}, {3, 4}, {3, 7}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{2, 4}, {3, 4}, {4, 8}}, 
     {{7, 8}, {3, 7}, {5, 7}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{7, 8}, {6, 8}, {2, 4}}, {{3, 4}, {7, 8}, {2, 4}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{2, 4}, {3, 4}, {4, 8}}}, 
    {{{2, 6}, {4, 8}, {3, 7}}, {{3, 7}, {1, 3}, {1, 2}}, 
     {{2, 6}, {3, 7}, {1, 2}}}, {{{3, 7}, {1, 3}, {1, 2}}, 
     {{7, 8}, {3, 7}, {1, 2}}, {{6, 8}, {7, 8}, {1, 2}}, 
     {{2, 6}, {6, 8}, {1, 2}}}, {{{1, 2}, {5, 7}, {1, 3}}, 
     {{1, 2}, {2, 6}, {5, 7}}, {{7, 8}, {5, 7}, {2, 6}}, 
     {{7, 8}, {2, 6}, {4, 8}}}, {{{1, 2}, {5, 7}, {1, 3}}, 
     {{6, 8}, {5, 7}, {2, 6}}, {{1, 2}, {2, 6}, {5, 7}}}, 
    {{{1, 2}, {5, 6}, {1, 3}}, {{6, 8}, {4, 8}, {3, 7}}, 
     {{6, 8}, {3, 7}, {5, 6}}, {{5, 6}, {3, 7}, {1, 3}}}, 
    {{{1, 2}, {5, 6}, {1, 3}}, {{7, 8}, {3, 7}, {5, 6}}, 
     {{5, 6}, {3, 7}, {1, 3}}}, {{{6, 8}, {4, 8}, {7, 8}}, 
     {{5, 7}, {1, 3}, {1, 2}}, {{5, 6}, {5, 7}, {1, 2}}}, 
    {{{5, 7}, {1, 3}, {1, 2}}, {{5, 6}, {5, 7}, {1, 2}}}, 
    {{{5, 6}, {4, 8}, {5, 7}}, {{5, 6}, {2, 6}, {4, 8}}, 
     {{5, 7}, {4, 8}, {3, 7}}, {{1, 2}, {1, 5}, {1, 3}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{5, 6}, {2, 6}, {6, 8}}, 
     {{7, 8}, {3, 7}, {5, 7}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{7, 8}, {5, 6}, {2, 6}}, {{4, 8}, {7, 8}, {2, 6}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{5, 6}, {2, 6}, {6, 8}}}, 
    {{{1, 2}, {1, 5}, {1, 3}}, {{4, 8}, {3, 7}, {5, 7}}, 
     {{6, 8}, {4, 8}, {5, 7}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{3, 7}, {5, 7}, {7, 8}}}, {{{1, 2}, {1, 5}, {1, 3}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{1, 2}, {1, 5}, {1, 3}}}, 
    {{{1, 5}, {1, 2}, {1, 3}}}, {{{1, 3}, {4, 8}, {1, 2}}, 
     {{1, 3}, {7, 8}, {4, 8}}, {{1, 3}, {1, 5}, {7, 8}}, 
     {{7, 8}, {1, 5}, {6, 8}}, {{1, 2}, {6, 8}, {1, 5}}, 
     {{1, 2}, {4, 8}, {6, 8}}}, {{{7, 8}, {5, 7}, {1, 5}}, 
     {{1, 2}, {7, 8}, {1, 5}}, {{3, 7}, {7, 8}, {1, 2}}, 
     {{1, 3}, {3, 7}, {1, 2}}}, {{{1, 3}, {4, 8}, {1, 2}}, 
     {{1, 3}, {3, 7}, {4, 8}}, {{6, 8}, {1, 2}, {4, 8}}, 
     {{6, 8}, {5, 7}, {1, 5}}, {{6, 8}, {1, 5}, {1, 2}}}, 
    {{{1, 3}, {1, 5}, {5, 6}}, {{6, 8}, {1, 3}, {5, 6}}, 
     {{1, 2}, {1, 3}, {6, 8}}, {{2, 6}, {1, 2}, {6, 8}}}, 
    {{{1, 2}, {1, 3}, {4, 8}}, {{1, 2}, {4, 8}, {2, 6}}, 
     {{5, 6}, {7, 8}, {1, 5}}, {{1, 3}, {1, 5}, {7, 8}}, 
     {{1, 3}, {7, 8}, {4, 8}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{6, 8}, {2, 6}, {1, 2}}, {{7, 8}, {6, 8}, {1, 2}}, 
     {{3, 7}, {7, 8}, {1, 2}}, {{1, 3}, {3, 7}, {1, 2}}}, 
    {{{1, 2}, {4, 8}, {2, 6}}, {{1, 2}, {1, 3}, {4, 8}}, 
     {{5, 6}, {5, 7}, {1, 5}}, {{1, 3}, {3, 7}, {4, 8}}}, 
    {{{5, 7}, {5, 6}, {1, 2}}, {{1, 3}, {5, 7}, {1, 2}}}, 
    {{{1, 2}, {1, 3}, {4, 8}}, {{1, 3}, {7, 8}, {4, 8}}, 
     {{1, 3}, {5, 7}, {7, 8}}, {{1, 2}, {6, 8}, {5, 6}}, 
     {{1, 2}, {4, 8}, {6, 8}}}, {{{1, 2}, {7, 8}, {5, 6}}, 
     {{1, 2}, {1, 3}, {3, 7}}, {{7, 8}, {1, 2}, {3, 7}}}, 
    {{{1, 2}, {6, 8}, {5, 6}}, {{1, 2}, {1, 3}, {4, 8}}, 
     {{1, 2}, {4, 8}, {6, 8}}, {{1, 3}, {3, 7}, {4, 8}}}, 
    {{{1, 3}, {5, 7}, {6, 8}}, {{6, 8}, {2, 6}, {1, 2}}, 
     {{1, 3}, {6, 8}, {1, 2}}}, {{{1, 3}, {5, 7}, {7, 8}}, 
     {{1, 2}, {4, 8}, {2, 6}}, {{1, 3}, {4, 8}, {1, 2}}, 
     {{1, 3}, {7, 8}, {4, 8}}}, {{{6, 8}, {2, 6}, {1, 2}}, 
     {{7, 8}, {6, 8}, {1, 2}}, {{3, 7}, {7, 8}, {1, 2}}, 
     {{1, 3}, {3, 7}, {1, 2}}}, {{{4, 8}, {1, 2}, {1, 3}}, 
     {{4, 8}, {1, 3}, {3, 7}}, {{4, 8}, {2, 6}, {1, 2}}}, 
    {{{4, 8}, {3, 4}, {1, 3}}, {{1, 5}, {4, 8}, {1, 3}}, 
     {{2, 4}, {4, 8}, {1, 5}}, {{1, 2}, {2, 4}, {1, 5}}}, 
    {{{1, 2}, {6, 8}, {1, 5}}, {{6, 8}, {7, 8}, {1, 5}}, 
     {{1, 2}, {2, 4}, {6, 8}}, {{3, 4}, {1, 3}, {7, 8}}, 
     {{7, 8}, {1, 3}, {1, 5}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{7, 8}, {5, 7}, {1, 5}}, {{4, 8}, {7, 8}, {1, 5}}, 
     {{2, 4}, {4, 8}, {1, 5}}, {{1, 2}, {2, 4}, {1, 5}}}, 
    {{{6, 8}, {5, 7}, {1, 5}}, {{3, 4}, {1, 3}, {3, 7}}, 
     {{6, 8}, {1, 2}, {2, 4}}, {{6, 8}, {1, 5}, {1, 2}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{4, 8}, {3, 4}, {1, 3}}, 
     {{6, 8}, {4, 8}, {1, 3}}, {{5, 6}, {6, 8}, {1, 3}}, 
     {{1, 5}, {5, 6}, {1, 3}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{5, 6}, {7, 8}, {1, 5}}, {{3, 4}, {1, 3}, {7, 8}}, 
     {{7, 8}, {1, 3}, {1, 5}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{5, 6}, {5, 7}, {1, 5}}, {{3, 4}, {1, 3}, {3, 7}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{5, 6}, {5, 7}, {1, 5}}, {{3, 4}, {1, 3}, {3, 7}}}, 
    {{{1, 2}, {2, 4}, {5, 6}}, {{5, 6}, {2, 4}, {4, 8}}, 
     {{3, 4}, {1, 3}, {5, 7}}, {{3, 4}, {5, 7}, {4, 8}}, 
     {{5, 6}, {4, 8}, {5, 7}}}, {{{6, 8}, {5, 6}, {1, 2}}, 
     {{2, 4}, {6, 8}, {1, 2}}, {{7, 8}, {3, 4}, {1, 3}}, 
     {{5, 7}, {7, 8}, {1, 3}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{1, 2}, {2, 4}, {5, 6}}, {{7, 8}, {5, 6}, {4, 8}}, 
     {{5, 6}, {2, 4}, {4, 8}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{5, 6}, {1, 2}, {2, 4}}, {{6, 8}, {5, 6}, {2, 4}}}, 
    {{{3, 4}, {1, 3}, {5, 7}}, {{1, 2}, {2, 4}, {2, 6}}, 
     {{6, 8}, {4, 8}, {5, 7}}, {{3, 4}, {5, 7}, {4, 8}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{7, 8}, {3, 4}, {1, 3}}, 
     {{5, 7}, {7, 8}, {1, 3}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{3, 4}, {1, 3}, {3, 7}}, {{7, 8}, {6, 8}, {4, 8}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{3, 4}, {1, 3}, {3, 7}}}, 
    {{{1, 5}, {1, 2}, {3, 4}}, {{3, 7}, {1, 5}, {3, 4}}}, 
    {{{7, 8}, {1, 5}, {6, 8}}, {{7, 8}, {3, 7}, {1, 5}}, 
     {{1, 2}, {6, 8}, {1, 5}}, {{1, 2}, {3, 4}, {4, 8}}, 
     {{1, 2}, {4, 8}, {6, 8}}}, {{{1, 2}, {3, 4}, {7, 8}}, 
     {{7, 8}, {5, 7}, {1, 5}}, {{1, 2}, {7, 8}, {1, 5}}}, 
    {{{6, 8}, {5, 7}, {1, 5}}, {{6, 8}, {1, 5}, {1, 2}}, 
     {{6, 8}, {1, 2}, {4, 8}}, {{1, 2}, {3, 4}, {4, 8}}}, 
    {{{1, 2}, {3, 4}, {2, 6}}, {{3, 4}, {6, 8}, {2, 6}}, 
     {{5, 6}, {3, 7}, {1, 5}}, {{5, 6}, {6, 8}, {3, 7}}, 
     {{3, 4}, {3, 7}, {6, 8}}}, {{{4, 8}, {2, 6}, {1, 2}}, 
     {{3, 4}, {4, 8}, {1, 2}}, {{3, 7}, {1, 5}, {5, 6}}, 
     {{7, 8}, {3, 7}, {5, 6}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{1, 2}, {3, 4}, {2, 6}}, {{6, 8}, {2, 6}, {3, 4}}, 
     {{6, 8}, {3, 4}, {7, 8}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{2, 6}, {1, 2}, {3, 4}}, {{4, 8}, {2, 6}, {3, 4}}}, 
    {{{1, 2}, {3, 4}, {5, 6}}, {{3, 4}, {3, 7}, {5, 7}}, 
     {{5, 6}, {3, 4}, {5, 7}}}, {{{1, 2}, {3, 4}, {4, 8}}, 
     {{1, 2}, {6, 8}, {5, 6}}, {{1, 2}, {4, 8}, {6, 8}}, 
     {{7, 8}, {3, 7}, {5, 7}}}, {{{7, 8}, {5, 6}, {1, 2}}, 
     {{3, 4}, {7, 8}, {1, 2}}}, {{{1, 2}, {3, 4}, {4, 8}}, 
     {{1, 2}, {6, 8}, {5, 6}}, {{1, 2}, {4, 8}, {6, 8}}}, 
    {{{6, 8}, {3, 7}, {5, 7}}, {{1, 2}, {3, 4}, {2, 6}}, 
     {{6, 8}, {2, 6}, {3, 4}}, {{6, 8}, {3, 4}, {3, 7}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}, {{4, 8}, {2, 6}, {1, 2}}, 
     {{3, 4}, {4, 8}, {1, 2}}}, {{{1, 2}, {3, 4}, {2, 6}}, 
     {{6, 8}, {2, 6}, {3, 4}}, {{6, 8}, {3, 4}, {7, 8}}}, 
    {{{4, 8}, {2, 6}, {1, 2}}, {{3, 4}, {4, 8}, {1, 2}}}, 
    {{{1, 5}, {4, 8}, {3, 7}}, {{1, 5}, {1, 2}, {2, 4}}, 
     {{4, 8}, {1, 5}, {2, 4}}}, {{{7, 8}, {3, 7}, {1, 5}}, 
     {{1, 2}, {2, 4}, {6, 8}}, {{1, 2}, {6, 8}, {1, 5}}, 
     {{7, 8}, {1, 5}, {6, 8}}}, {{{1, 5}, {1, 2}, {2, 4}}, 
     {{5, 7}, {1, 5}, {2, 4}}, {{7, 8}, {5, 7}, {2, 4}}, 
     {{4, 8}, {7, 8}, {2, 4}}}, {{{6, 8}, {1, 2}, {2, 4}}, 
     {{6, 8}, {5, 7}, {1, 5}}, {{6, 8}, {1, 5}, {1, 2}}}, 
    {{{5, 6}, {6, 8}, {3, 7}}, {{5, 6}, {3, 7}, {1, 5}}, 
     {{6, 8}, {4, 8}, {3, 7}}, {{1, 2}, {2, 4}, {2, 6}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{3, 7}, {1, 5}, {5, 6}}, 
     {{7, 8}, {3, 7}, {5, 6}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{5, 6}, {5, 7}, {1, 5}}, {{7, 8}, {6, 8}, {4, 8}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{5, 6}, {5, 7}, {1, 5}}}, 
    {{{5, 7}, {5, 6}, {4, 8}}, {{5, 7}, {4, 8}, {3, 7}}, 
     {{1, 2}, {2, 4}, {5, 6}}, {{5, 6}, {2, 4}, {4, 8}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}, {{5, 6}, {1, 2}, {2, 4}}, 
     {{6, 8}, {5, 6}, {2, 4}}}, {{{1, 2}, {2, 4}, {5, 6}}, 
     {{7, 8}, {5, 6}, {4, 8}}, {{5, 6}, {2, 4}, {4, 8}}}, 
    {{{5, 6}, {1, 2}, {2, 4}}, {{6, 8}, {5, 6}, {2, 4}}}, 
    {{{1, 2}, {2, 4}, {2, 6}}, {{4, 8}, {3, 7}, {5, 7}}, 
     {{6, 8}, {4, 8}, {5, 7}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{7, 8}, {3, 7}, {5, 7}}}, {{{1, 2}, {2, 4}, {2, 6}}, 
     {{7, 8}, {6, 8}, {4, 8}}}, {{{1, 2}, {2, 4}, {2, 6}}}, 
    {{{2, 6}, {2, 4}, {1, 3}}, {{1, 5}, {2, 6}, {1, 3}}}, 
    {{{2, 4}, {1, 3}, {4, 8}}, {{1, 3}, {7, 8}, {4, 8}}, 
     {{6, 8}, {1, 5}, {2, 6}}, {{6, 8}, {7, 8}, {1, 5}}, 
     {{1, 3}, {1, 5}, {7, 8}}}, {{{2, 4}, {7, 8}, {2, 6}}, 
     {{5, 7}, {2, 6}, {7, 8}}, {{5, 7}, {1, 5}, {2, 6}}, 
     {{2, 4}, {1, 3}, {3, 7}}, {{2, 4}, {3, 7}, {7, 8}}}, 
    {{{4, 8}, {2, 4}, {1, 3}}, {{3, 7}, {4, 8}, {1, 3}}, 
     {{1, 5}, {2, 6}, {6, 8}}, {{5, 7}, {1, 5}, {6, 8}}}, 
    {{{2, 4}, {1, 3}, {6, 8}}, {{1, 3}, {1, 5}, {5, 6}}, 
     {{6, 8}, {1, 3}, {5, 6}}}, {{{5, 6}, {7, 8}, {1, 5}}, 
     {{2, 4}, {1, 3}, {4, 8}}, {{7, 8}, {4, 8}, {1, 3}}, 
     {{7, 8}, {1, 3}, {1, 5}}}, {{{2, 4}, {7, 8}, {6, 8}}, 
     {{2, 4}, {1, 3}, {3, 7}}, {{2, 4}, {3, 7}, {7, 8}}, 
     {{5, 6}, {5, 7}, {1, 5}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{4, 8}, {2, 4}, {1, 3}}, {{3, 7}, {4, 8}, {1, 3}}}, 
    {{{2, 4}, {1, 3}, {5, 7}}, {{5, 7}, {5, 6}, {2, 6}}, 
     {{2, 4}, {5, 7}, {2, 6}}}, {{{1, 3}, {7, 8}, {4, 8}}, 
     {{1, 3}, {5, 7}, {7, 8}}, {{1, 3}, {4, 8}, {2, 4}}, 
     {{5, 6}, {2, 6}, {6, 8}}}, {{{2, 4}, {1, 3}, {3, 7}}, 
     {{2, 4}, {3, 7}, {7, 8}}, {{5, 6}, {2, 6}, {7, 8}}, 
     {{2, 4}, {7, 8}, {2, 6}}}, {{{5, 6}, {2, 6}, {6, 8}}, 
     {{4, 8}, {2, 4}, {1, 3}}, {{3, 7}, {4, 8}, {1, 3}}}, 
    {{{6, 8}, {2, 4}, {1, 3}}, {{5, 7}, {6, 8}, {1, 3}}}, 
    {{{1, 3}, {7, 8}, {4, 8}}, {{1, 3}, {5, 7}, {7, 8}}, 
     {{1, 3}, {4, 8}, {2, 4}}}, {{{2, 4}, {7, 8}, {6, 8}}, 
     {{2, 4}, {1, 3}, {3, 7}}, {{2, 4}, {3, 7}, {7, 8}}}, 
    {{{4, 8}, {2, 4}, {1, 3}}, {{3, 7}, {4, 8}, {1, 3}}}, 
    {{{2, 6}, {4, 8}, {1, 5}}, {{4, 8}, {3, 4}, {1, 3}}, 
     {{1, 5}, {4, 8}, {1, 3}}}, {{{3, 4}, {1, 3}, {7, 8}}, 
     {{6, 8}, {1, 5}, {2, 6}}, {{6, 8}, {7, 8}, {1, 5}}, 
     {{7, 8}, {1, 3}, {1, 5}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{7, 8}, {2, 6}, {4, 8}}, {{5, 7}, {2, 6}, {7, 8}}, 
     {{5, 7}, {1, 5}, {2, 6}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{1, 5}, {2, 6}, {6, 8}}, {{5, 7}, {1, 5}, {6, 8}}}, 
    {{{4, 8}, {3, 4}, {1, 3}}, {{6, 8}, {4, 8}, {1, 3}}, 
     {{5, 6}, {6, 8}, {1, 3}}, {{1, 5}, {5, 6}, {1, 3}}}, 
    {{{7, 8}, {1, 5}, {5, 6}}, {{7, 8}, {3, 4}, {1, 3}}, 
     {{7, 8}, {1, 3}, {1, 5}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{3, 4}, {1, 3}, {3, 7}}, {{6, 8}, {4, 8}, {7, 8}}}, 
    {{{5, 6}, {5, 7}, {1, 5}}, {{3, 4}, {1, 3}, {3, 7}}}, 
    {{{3, 4}, {1, 3}, {5, 7}}, {{5, 6}, {2, 6}, {4, 8}}, 
     {{5, 6}, {4, 8}, {5, 7}}, {{3, 4}, {5, 7}, {4, 8}}}, 
    {{{5, 6}, {2, 6}, {6, 8}}, {{7, 8}, {3, 4}, {1, 3}}, 
     {{5, 7}, {7, 8}, {1, 3}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{7, 8}, {5, 6}, {2, 6}}, {{4, 8}, {7, 8}, {2, 6}}}, 
    {{{3, 4}, {1, 3}, {3, 7}}, {{5, 6}, {2, 6}, {6, 8}}}, 
    {{{3, 4}, {1, 3}, {5, 7}}, {{6, 8}, {4, 8}, {5, 7}}, 
     {{3, 4}, {5, 7}, {4, 8}}}, {{{7, 8}, {3, 4}, {1, 3}}, 
     {{5, 7}, {7, 8}, {1, 3}}}, {{{3, 4}, {1, 3}, {3, 7}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{3, 4}, {1, 3}, {3, 7}}}, 
    {{{2, 6}, {3, 7}, {1, 5}}, {{2, 6}, {2, 4}, {3, 4}}, 
     {{3, 7}, {2, 6}, {3, 4}}}, {{{7, 8}, {1, 5}, {6, 8}}, 
     {{7, 8}, {3, 7}, {1, 5}}, {{6, 8}, {1, 5}, {2, 6}}, 
     {{2, 4}, {3, 4}, {4, 8}}}, {{{5, 7}, {1, 5}, {2, 6}}, 
     {{2, 4}, {3, 4}, {7, 8}}, {{2, 4}, {7, 8}, {2, 6}}, 
     {{5, 7}, {2, 6}, {7, 8}}}, {{{2, 4}, {3, 4}, {4, 8}}, 
     {{1, 5}, {2, 6}, {6, 8}}, {{5, 7}, {1, 5}, {6, 8}}}, 
    {{{5, 6}, {3, 7}, {1, 5}}, {{2, 4}, {3, 4}, {6, 8}}, 
     {{3, 7}, {6, 8}, {3, 4}}, {{5, 6}, {6, 8}, {3, 7}}}, 
    {{{2, 4}, {3, 4}, {4, 8}}, {{3, 7}, {1, 5}, {5, 6}}, 
     {{7, 8}, {3, 7}, {5, 6}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{6, 8}, {2, 4}, {3, 4}}, {{7, 8}, {6, 8}, {3, 4}}}, 
    {{{5, 6}, {5, 7}, {1, 5}}, {{2, 4}, {3, 4}, {4, 8}}}, 
    {{{3, 4}, {3, 7}, {5, 7}}, {{2, 4}, {3, 4}, {5, 7}}, 
     {{2, 6}, {2, 4}, {5, 7}}, {{5, 6}, {2, 6}, {5, 7}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}, {{5, 6}, {2, 6}, {6, 8}}, 
     {{2, 4}, {3, 4}, {4, 8}}}, {{{7, 8}, {2, 4}, {3, 4}}, 
     {{7, 8}, {5, 6}, {2, 6}}, {{7, 8}, {2, 6}, {2, 4}}}, 
    {{{5, 6}, {2, 6}, {6, 8}}, {{2, 4}, {3, 4}, {4, 8}}}, 
    {{{6, 8}, {3, 7}, {5, 7}}, {{6, 8}, {2, 4}, {3, 4}}, 
     {{6, 8}, {3, 4}, {3, 7}}}, {{{7, 8}, {3, 7}, {5, 7}}, 
     {{2, 4}, {3, 4}, {4, 8}}}, {{{7, 8}, {6, 8}, {2, 4}}, 
     {{3, 4}, {7, 8}, {2, 4}}}, {{{2, 4}, {3, 4}, {4, 8}}}, 
    {{{3, 7}, {1, 5}, {2, 6}}, {{4, 8}, {3, 7}, {2, 6}}}, 
    {{{7, 8}, {3, 7}, {1, 5}}, {{7, 8}, {1, 5}, {6, 8}}, 
     {{1, 5}, {2, 6}, {6, 8}}}, {{{5, 7}, {1, 5}, {2, 6}}, 
     {{5, 7}, {2, 6}, {7, 8}}, {{2, 6}, {4, 8}, {7, 8}}}, 
    {{{1, 5}, {2, 6}, {6, 8}}, {{5, 7}, {1, 5}, {6, 8}}}, 
    {{{5, 6}, {3, 7}, {1, 5}}, {{5, 6}, {6, 8}, {3, 7}}, 
     {{6, 8}, {4, 8}, {3, 7}}}, {{{3, 7}, {1, 5}, {5, 6}}, 
     {{7, 8}, {3, 7}, {5, 6}}}, {{{5, 6}, {5, 7}, {1, 5}}, 
     {{6, 8}, {4, 8}, {7, 8}}}, {{{5, 6}, {5, 7}, {1, 5}}}, 
    {{{4, 8}, {5, 7}, {5, 6}}, {{4, 8}, {3, 7}, {5, 7}}, 
     {{4, 8}, {5, 6}, {2, 6}}}, {{{7, 8}, {3, 7}, {5, 7}}, 
     {{5, 6}, {2, 6}, {6, 8}}}, {{{4, 8}, {7, 8}, {5, 6}}, 
     {{2, 6}, {4, 8}, {5, 6}}}, {{{5, 6}, {2, 6}, {6, 8}}}, 
    {{{3, 7}, {5, 7}, {6, 8}}, {{4, 8}, {3, 7}, {6, 8}}}, 
    {{{7, 8}, {3, 7}, {5, 7}}}, {{{6, 8}, {4, 8}, {7, 8}}}, {}}; 

    std::vector<std::vector<std::vector<std::pair<int,int>>>> pyramidCycles(32);

    pyramidCycles = {{}, {{{1, 5}, {3, 5}, {4, 5}}}, {{{1, 4}, {2, 4}, {4, 5}}}, {{{1,
    4}, {1, 5}, {3, 5}, {2, 4}}}, {{{1, 3}, {2, 3}, {3, 5}}}, {{{1,
    3}, {2, 3}, {4, 5}, {1, 5}}}, {{{3, 5}, {1, 3}, {1, 4}, {4,
    5}}, {{1, 3}, {2, 3}, {2, 4}, {1, 4}}}, {{{1, 3}, {2, 3}, {2,
    4}, {1, 4}, {1, 5}}}, {{{1, 2}, {2, 3}, {2, 4}}}, {{{1, 2}, {2,
    4}, {4, 5}, {1, 5}}, {{1, 2}, {2, 3}, {3, 5}, {1, 5}}}, {{{1,
    2}, {2, 3}, {4, 5}, {1, 4}}}, {{{1, 2}, {1, 4}, {1, 5}, {3,
    5}, {2, 3}}}, {{{1, 2}, {2, 4}, {3, 5}, {1, 3}}}, {{{1, 2}, {1,
    3}, {1, 5}, {4, 5}, {2, 4}}}, {{{1, 4}, {1, 2}, {1, 3}, {3,
    5}, {4, 5}}}, {{{1, 2}, {1, 3}, {1, 5}, {1, 4}}}, {{{1, 3}, {1,
    2}, {1, 4}, {1, 5}}}, {{{1, 4}, {1, 2}, {1, 3}, {3, 5}, {4,
    5}}}, {{{1, 2}, {1, 3}, {1, 5}, {4, 5}, {2, 4}}}, {{{1, 2}, {1,
    3}, {3, 5}, {2, 4}}}, {{{1, 2}, {1, 4}, {1, 5}, {3, 5}, {2,
    3}}}, {{{1, 2}, {2, 3}, {4, 5}, {1, 4}}}, {{{1, 2}, {2, 3}, {2,
    4}}, {{1, 5}, {4, 5}, {3, 5}}}, {{{1, 2}, {2, 3}, {2, 4}}}, {{{1,
    3}, {1, 5}, {1, 4}, {2, 4}, {2, 3}}}, {{{1, 3}, {2, 3}, {3,
    5}}, {{1, 4}, {2, 4}, {4, 5}}}, {{{1, 3}, {2, 3}, {4, 5}, {1,
    5}}}, {{{1, 3}, {2, 3}, {3, 5}}}, {{{1, 4}, {1, 5}, {3, 5}, {2,
    4}}}, {{{1, 4}, {2, 4}, {4, 5}}}, {{{1, 5}, {3, 5}, {4, 5}}}, {}};

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
    
    auto getLinearRoot = [](const Point3D& v1, const Point3D& v2, double val1, double val2) -> Point3D {
        // Avoid division by zero - if values are same, return midpoint
        if (std::abs(val1 - val2) < 1e-12) {
            return Point3D((v1.x + v2.x) * 0.5, (v1.y + v2.y) * 0.5, (v1.z + v2.z) * 0.5);
        }
        
        // Linear interpolation: t such that val1 + t*(val2 - val1) = 0
        double t = -val1 / (val2 - val1);
        
        // Clamp t to [0.05, 0.95] to leave 0.05 padding on each side
        //t = std::max(0.1, std::min(0.9, t));
        
        return Point3D(
            v1.x + t * (v2.x - v1.x),
            v1.y + t * (v2.y - v1.y),
            v1.z + t * (v2.z - v1.z)
        );
    };

    int zeroVertCounter = 0;
    for (const auto& vIdxs : cells) {

        const int n = static_cast<int>(vIdxs.size());
        if (n!=4 && n!=5 && n!=8){
            std::cout << "Degenerate cell with " << n << " vertices" << std::endl;
            continue;
        }  // degenerate


        std::vector<int>signs(n);
        for(int i = 0; i < n; ++i){
            signs[i] = getSign(vals[vIdxs[i]]);
        }

        if (n == 4){
            
            int patternIdx = signs[0] * 8 + signs[1] * 4 + signs[2] * 2 + signs[3] * 1;
            
            // Get lookup table entry
            const auto& entry = tetTableEinds[patternIdx];
            
            // Skip if no zero-crossing (empty entry)
            if (entry.edgeIndices.empty()) {
                continue;
            }

            int ori = getTetOrientation(
                vts[vIdxs[0]], vts[vIdxs[1]], 
                vts[vIdxs[2]], vts[vIdxs[3]]
            );

            // Create zero-crossing points for each edge in the entry
            // Store local edge index -> zero vertex index mapping
            std::vector<int> localEdgeToZeroVert(6, -1);

            for (int edgeIdx1Based : entry.edgeIndices) {
                // Convert from 1-based to 0-based edge index
                int edgeIdx = edgeIdx1Based - 1;
                
                // Get vertex pair for this edge
                const auto& edgePair = tetEdges[edgeIdx];
                int vIdx1 = vIdxs[edgePair.first];
                int vIdx2 = vIdxs[edgePair.second];
                
                // Verify that signs differ (should always be true for edges in the lookup table)
                int sign1 = getSign(vals[vIdx1]);
                int sign2 = getSign(vals[vIdx2]);
                if (sign1 == sign2) {
                    std::cout<<"THIS HAPPENED" <<std::endl;
                    // This shouldn't happen if lookup table is correct, but skip to be safe
                    continue;
                }
                
                // Create sorted edge key for hashing
                UndirectedEdgeKey edgeKey(vIdx1, vIdx2);
                
                // Check if we already have a zero point for this edge
                auto it = edgeHash.find(edgeKey);
                if (it != edgeHash.end()) {
                    localEdgeToZeroVert[edgeIdx] = it->second;
                } else {

                    Point3D zeroPt;

                    if (vIdx1 >= originalNumVertices || vIdx2 >= originalNumVertices || (vIdx1< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx1] )|| (vIdx2< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx2])) {

                        zeroPt.x = 0.5 * (vts[vIdx1].x + vts[vIdx2].x);
                        zeroPt.y = 0.5 * (vts[vIdx1].y + vts[vIdx2].y);
                        zeroPt.z = 0.5 * (vts[vIdx1].z + vts[vIdx2].z);
                        usedMidPoint.push_back(true);
                    }
                    else{
                        zeroPt = getLinearRoot(
                            vts[vIdx1], vts[vIdx2],
                            vals[vIdx1], vals[vIdx2]
                        );
                        usedMidPoint.push_back(false);
                    }
                    // Compute zero-crossing point

                    // Add to zero vertices
                    int zeroVertIdx = zeroVertCounter++;
                    zeroVerts.push_back(zeroPt);
                    edgeHash[edgeKey] = zeroVertIdx;
                    localEdgeToZeroVert[edgeIdx] = zeroVertIdx;
                    zeroCrossingEdges.push_back(std::make_pair(vts[vIdx1], vts[vIdx2]));
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
                        // Only add if not already in the list (avoid duplicates)
                        if (std::find(triVerts.begin(), triVerts.end(), zeroVertIdx) == triVerts.end()) {
                            triVerts.push_back(zeroVertIdx);
                        }
                    }
                }
                
                // Only create triangle if we have exactly 3 distinct vertices
                if (triVerts.size() == 3) {
                    // Check that all three vertices are distinct
                    if (triVerts[0] != triVerts[1] && triVerts[1] != triVerts[2] && triVerts[0] != triVerts[2]) {
                        if (ori == 1) {
                            // Positive orientation: use cycle order as-is
                            zeroTris.emplace_back(triVerts[0], triVerts[1], triVerts[2]);

                        } else {
                            // Negative or zero orientation: reverse the cycle order
                            zeroTris.emplace_back(triVerts[2], triVerts[1], triVerts[0]);
                        }
                    } else {
                        std::cout << "Skipping degenerate triangle with duplicate vertices: [" 
                                  << triVerts[0] << ", " << triVerts[1] << ", " << triVerts[2] << "]" << std::endl;
                    }
                }
                else{
                    std::cout << "Degenerate triangle with " << triVerts.size() << " vertices" << std::endl;
                }
            }
        }
        else if (n == 5){
            
            int key = signs[0] * 16 + signs[1] * 8 + signs[2] * 4 + signs[3] * 2 + signs[4] * 1;

            if (key < pyramidCycles.size() && !pyramidCycles[key].empty()) {
                for (int tri_idx = 0; tri_idx < pyramidCycles[key].size(); ++tri_idx) {
                    auto& triangle = pyramidCycles[key][tri_idx];
            
                    std::vector<int> triangleVertices;
                    triangleVertices.reserve(5);
                    
                    for (const auto& edge : triangle) {
                        // Convert from 1-based to 0-based edge index
                        int vIdx1 = vIdxs[edge.first-1];
                        int vIdx2 = vIdxs[edge.second-1];

                        int sign1 = getSign(vals[vIdx1]);
                        int sign2 = getSign(vals[vIdx2]);
                        if (sign1 == sign2) {
                            // This shouldn't happen if lookup table is correct, but skip to be safe
                            continue;
                        }
                        
                        // Create sorted edge key for hashing
                        UndirectedEdgeKey edgeKey(vIdx1, vIdx2);
                        
                        // Check if we already have a zero point for this edge
                        auto it = edgeHash.find(edgeKey);
                        int zeroVertIdx;
                        if (it != edgeHash.end()) {
                            zeroVertIdx = it->second;
                        } else {

                            Point3D zeroPt;

                            if (vIdx1 >= originalNumVertices || vIdx2 >= originalNumVertices || (vIdx1< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx1] )|| (vIdx2< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx2])) {

                                zeroPt.x = 0.5 * (vts[vIdx1].x + vts[vIdx2].x);
                                zeroPt.y = 0.5 * (vts[vIdx1].y + vts[vIdx2].y);
                                zeroPt.z = 0.5 * (vts[vIdx1].z + vts[vIdx2].z);
                                usedMidPoint.push_back(true);
                            }
                            else{
                                // Compute zero-crossing point
                                zeroPt = getLinearRoot(
                                    vts[vIdx1], vts[vIdx2],
                                    vals[vIdx1], vals[vIdx2]
                                );
                                usedMidPoint.push_back(false);
                            }
                            
                            // Add to zero vertices
                            zeroVertIdx = zeroVertCounter++;
                            zeroVerts.push_back(zeroPt);
                            edgeHash[edgeKey] = zeroVertIdx;
                            zeroCrossingEdges.push_back(std::make_pair(vts[vIdx1], vts[vIdx2]));
                        }
                        // Only add if not already in the list (avoid duplicates)
                        if (std::find(triangleVertices.begin(), triangleVertices.end(), zeroVertIdx) == triangleVertices.end()) {
                            triangleVertices.push_back(zeroVertIdx);
                        }
                    }
                    if (triangleVertices.size() == 3) {
                        // Check that all three vertices are distinct
                        if (triangleVertices[0] != triangleVertices[1] && 
                            triangleVertices[1] != triangleVertices[2] && 
                            triangleVertices[0] != triangleVertices[2]) {
                            zeroTris.emplace_back(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
                        } else {
                            std::cout << "Skipping degenerate triangle (n==5) with duplicate vertices: [" 
                                      << triangleVertices[0] << ", " << triangleVertices[1] << ", " << triangleVertices[2] << "]" << std::endl;
                        }
                    }
                    else if (triangleVertices.size()>=3){
                        //this is a fan. use first vertex as pivot and add the rest of the vertices to the fan
                        int pivotIdx = triangleVertices[0];
                        for (int i = 1; i < triangleVertices.size(); ++i) {
                            int v1 = pivotIdx;
                            int v2 = triangleVertices[i];
                            int v3 = triangleVertices[(i+1)%triangleVertices.size()];
                            // Only add triangle if all three vertices are distinct
                            if (v1 != v2 && v2 != v3 && v1 != v3) {
                                zeroTris.emplace_back(v1, v2, v3);
                            }
                        }
                    }
                    
                }
            }
        }
        else if (n == 8){
            
            int key = signs[0] * 128 + signs[1] * 64 + signs[2] * 32 + signs[3] * 16+ signs[4] * 8 + signs[5] * 4 + signs[6] * 2 + signs[7] * 1;

            if (key < cubeCycles.size() && !cubeCycles[key].empty()) {
                for (int tri_idx = 0; tri_idx < cubeCycles[key].size(); ++tri_idx) {
                    auto& triangle = cubeCycles[key][tri_idx];

                    std::vector<int> triangleVertices;
                    triangleVertices.reserve(3);

                    for (const auto& edge : triangle) {
                        // Convert from 1-based to 0-based edge index
                        int vIdx1 = vIdxs[edge.first-1];
                        int vIdx2 = vIdxs[edge.second-1];

                        int sign1 = getSign(vals[vIdx1]);
                        int sign2 = getSign(vals[vIdx2]);
                        if (sign1 == sign2) {
                            // This shouldn't happen if lookup table is correct, but skip to be safe
                            continue;
                        }
                        
                        // Create sorted edge key for hashing
                        UndirectedEdgeKey edgeKey(vIdx1, vIdx2);
                        
                        // Check if we already have a zero point for this edge
                        auto it = edgeHash.find(edgeKey);
                        int zeroVertIdx;
                        if (it != edgeHash.end()) {
                            zeroVertIdx = it->second;
                        } else {
                            
                            Point3D zeroPt;

                            if (vIdx1 >= originalNumVertices || vIdx2 >= originalNumVertices || (vIdx1< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx1] )|| (vIdx2< vertNeedMidPoint.size()&&vertNeedMidPoint[vIdx2])) {

                                zeroPt.x = 0.5 * (vts[vIdx1].x + vts[vIdx2].x);
                                zeroPt.y = 0.5 * (vts[vIdx1].y + vts[vIdx2].y);
                                zeroPt.z = 0.5 * (vts[vIdx1].z + vts[vIdx2].z);
                                usedMidPoint.push_back(true);
                            }
                            else{
                                zeroPt = getLinearRoot(
                                    vts[vIdx1], vts[vIdx2],
                                    vals[vIdx1], vals[vIdx2]
                                );
                                usedMidPoint.push_back(false);
                            }
                            
                            // Add to zero vertices
                            zeroVertIdx = zeroVertCounter++;
                            zeroVerts.push_back(zeroPt);
                            edgeHash[edgeKey] = zeroVertIdx;
                            zeroCrossingEdges.push_back(std::make_pair(vts[vIdx1], vts[vIdx2]));
                        }
                        // Only add if not already in the list (avoid duplicates)
                        if (std::find(triangleVertices.begin(), triangleVertices.end(), zeroVertIdx) == triangleVertices.end()) {
                            triangleVertices.push_back(zeroVertIdx);
                        }
                    }
                    if (triangleVertices.size() == 3) {
                        // Check that all three vertices are distinct
                        if (triangleVertices[0] != triangleVertices[1] && 
                            triangleVertices[1] != triangleVertices[2] && 
                            triangleVertices[0] != triangleVertices[2]) {
                            zeroTris.emplace_back(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
                        } else {
                            std::cout << "Skipping degenerate triangle (n==8) with duplicate vertices: [" 
                                      << triangleVertices[0] << ", " << triangleVertices[1] << ", " << triangleVertices[2] << "]" << std::endl;
                        }
                    }
                }

            }
        }
    } // for cells

    // Clear edgeHash to free memory (no longer needed)
    edgeHash.clear();
    edgeHash.rehash(0);  // Shrink hash table
    
    // Remove duplicate triangles
    // Normalize triangles by sorting vertices, then use a set to track unique ones
    struct TriangleKey {
        int v1, v2, v3;
        TriangleKey(int a, int b, int c) {
            // Sort vertices to normalize triangle representation
            if (a > b) std::swap(a, b);
            if (b > c) std::swap(b, c);
            if (a > b) std::swap(a, b);
            v1 = a; v2 = b; v3 = c;
        }
        bool operator==(const TriangleKey& other) const {
            return v1 == other.v1 && v2 == other.v2 && v3 == other.v3;
        }
    };
    
    struct TriangleKeyHash {
        std::size_t operator()(const TriangleKey& t) const {
            return std::hash<int>()(t.v1) ^ (std::hash<int>()(t.v2) << 1) ^ (std::hash<int>()(t.v3) << 2);
        }
    };
    
    std::unordered_set<TriangleKey, TriangleKeyHash> uniqueTriangles;
    std::vector<Triangle> deduplicatedTris;
    deduplicatedTris.reserve(zeroTris.size());
    
    int duplicateCount = 0;
    for (const auto& tri : zeroTris) {
        TriangleKey key(tri.v1, tri.v2, tri.v3);
        if (uniqueTriangles.find(key) == uniqueTriangles.end()) {
            uniqueTriangles.insert(key);
            deduplicatedTris.push_back(tri);
        } else {
            duplicateCount++;
        }
    }
    
    if (duplicateCount > 0) {
        std::cout << "Removed " << duplicateCount << " duplicate triangles" << std::endl;
    }
    
    zeroTris = std::move(deduplicatedTris);
    
    // Shrink vectors to actual size
    zeroVerts.shrink_to_fit();
    zeroTris.shrink_to_fit();
    
    std::cout << "Memory optimized: cleared edgeHash, shrunk vectors" << std::endl;

    std::vector<std::vector<int>> zeroTrisVerts;
    for (const auto& t : zeroTris) {
        zeroTrisVerts.push_back({t.v1, t.v2, t.v3});
    }
    /*do fairing here*/
    {
        // Step 1: Build incidentVerts - for each vertex, record indices of vertices that share a triangle with it
        std::vector<std::vector<int>> incidentVerts(zeroVerts.size());
        for (const auto& tri : zeroTrisVerts) {
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
        std::vector<std::pair<Point3D, Point3D>> transformedEdges = zeroCrossingEdges;
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
        
        // Step 3: Iterative fairing (10 iterations)
        const int numIterations = 10;
        for (int iter = 0; iter < numIterations; ++iter) {
            std::vector<Point3D> newVerts = zeroVerts; // Temporary storage
            
            // Process each vertex that used midpoint
            for (size_t i = 0; i < zeroVerts.size(); ++i) {
                if (usedMidPoint[i]) {
                    // Compute 3D centroid of incident vertices
                    if (incidentVerts[i].empty()) {
                        continue; // Skip if no incident vertices
                    }
                    
                    Point3D centroid(0, 0, 0);
                    for (int j : incidentVerts[i]) {
                        centroid.x += zeroVerts[j].x;
                        centroid.y += zeroVerts[j].y;
                        centroid.z += zeroVerts[j].z;
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
            zeroVerts = newVerts;
        }
        
        std::cout << "Fairing complete: " << numIterations << " iterations applied" << std::endl;
    }

    
    std::ofstream objFile("output/output_surface.obj");
    if (objFile.is_open()) {
        for (const auto& vertex : zeroVerts) {
            objFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
        }
        
        int largeTriangleCount = 0;
        int invalidTriangleCount = 0;
        
        for (size_t triIdx = 0; triIdx < zeroTris.size(); ++triIdx) {
            const auto& triangle = zeroTris[triIdx];
            
            // OBJ format uses 1-based indexing for faces
            objFile << "f " << (triangle.v1 + 1) << " " << (triangle.v2 + 1) << " " << (triangle.v3 + 1) << std::endl;
        }
        
        if (largeTriangleCount > 0) {
            std::cout << "Found " << largeTriangleCount << " triangles with edge length > 5.0" << std::endl;
        }
        if (invalidTriangleCount > 0) {
            std::cerr << "Found " << invalidTriangleCount << " triangles with invalid vertex indices" << std::endl;
        }
    }
    objFile.close();
    std::cout << "Wrote " << zeroVerts.size() << " vertices and " << zeroTris.size()
              << " triangles to output_surface.obj in .obj format" << std::endl;


    return {zeroVerts, zeroTris};
}

// Function to create a higher resolution 3D array (that takes both allSimplices and original data)
std::vector<std::vector<std::vector<double>>> createHighResArray(
    const std::vector<Simplex*>& allSimplices,
    std::vector<std::vector<std::vector<double>>>& data,
    int dims[3]) {

    // Create higher resolution array (2x in each dimension)
    int highResDims[3] = {dims[0] * 2 - 1, dims[1] * 2 - 1, dims[2] * 2 - 1};
    
    std::vector<std::vector<std::vector<double>>> newCCForMarked(highResDims[0]);

    for (int i = 0; i < highResDims[0]; ++i) {
        newCCForMarked[i].resize(highResDims[1]);
        for (int j = 0; j < highResDims[1]; ++j) {
            newCCForMarked[i][j].resize(highResDims[2], -10000);
        }
    }

    std::cout << "Finish initing place" << std::endl;

    for (int z = 0; z < highResDims[2]; ++z) {
        for (int y = 0; y < highResDims[1]; ++y) {
            for (int x = 0; x < highResDims[0]; ++x) {

                int xOrig = x / 2;
                int yOrig = y / 2;
                int zOrig = z / 2;

                int xMod = x % 2;
                int yMod = y % 2;
                int zMod = z % 2;

                if (xMod == 0 && yMod == 0 && zMod == 0) {
        
                    newCCForMarked[x][y][z] = data[xOrig][yOrig][zOrig];
                }
                
                else if (xMod == 0 && yMod == 0 && zMod == 1) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig][yOrig][zOrig+1]});
                }
                else if (xMod == 0 && yMod == 1 && zMod == 0) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig][yOrig+1][zOrig]});
                }
                else if (xMod == 1 && yMod == 0 && zMod == 0) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig+1][yOrig][zOrig]});
                }
                else if (xMod == 1 && yMod == 0 && zMod == 1) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig+1][yOrig][zOrig], data[xOrig][yOrig][zOrig+1], data[xOrig+1][yOrig][zOrig+1]});
                }
                
                else if (xMod == 1 && yMod == 1 && zMod == 0) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig][yOrig+1][zOrig], data[xOrig+1][yOrig][zOrig], data[xOrig+1][yOrig+1][zOrig]});
                }
                else if (xMod == 0 && yMod == 1 && zMod == 1) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig][yOrig+1][zOrig], data[xOrig][yOrig][zOrig+1], data[xOrig][yOrig+1][zOrig+1]});
                }
                
                else if (xMod == 1 && yMod == 1 && zMod == 1) {
                    newCCForMarked[x][y][z] = std::max({data[xOrig][yOrig][zOrig], data[xOrig][yOrig+1][zOrig], data[xOrig+1][yOrig][zOrig], data[xOrig+1][yOrig+1][zOrig], data[xOrig][yOrig][zOrig+1], data[xOrig][yOrig+1][zOrig+1], data[xOrig+1][yOrig][zOrig+1], data[xOrig+1][yOrig+1][zOrig+1]});
                }
                

            }
        }
    }

    std::cout << "Finish populating newCCForMarked" << std::endl;

    //the size - 1 already removes the infinite cube
    for (int i = 0; i < allSimplices.size(); ++i) {
        Simplex* simplex = allSimplices[i];
        //if its val is infinity, skip it
        if (simplex->val == std::numeric_limits<double>::infinity()) {
            std::cout << "Skipping infinite Simplice " << i << std::endl;
            continue;
        }
        double alpha = simplex->val;
        
        int x = simplex->x;
        int y = simplex->y;
        int z = simplex->z;

        newCCForMarked[x][y][z] = alpha;
        
    }

    return newCCForMarked;
}


// Compute marked cubes from the high-resolution field using edge/face/cube sign-change criteria
std::vector<std::vector<std::vector<bool>>> computeMarkedFromHighRes(
    const std::vector<std::vector<std::vector<double>>>& newCCForMarked) {
    const int hx = static_cast<int>(newCCForMarked.size());
    const int hy = static_cast<int>(newCCForMarked[0].size());
    const int hz = static_cast<int>(newCCForMarked[0][0].size());

    // original dims from high-res dims: highRes = 2*orig - 1 => orig = (highRes + 1)/2
    const int ox = (hx + 1) / 2;
    const int oy = (hy + 1) / 2;
    const int oz = (hz + 1) / 2;

    // marked dims = original dims - 1
    std::vector<std::vector<std::vector<bool>>> marked(ox - 1);
    for (int i = 0; i < ox - 1; ++i) {
        marked[i].resize(oy - 1);
        for (int j = 0; j < oy - 1; ++j) {
            marked[i][j].resize(oz - 1, false);
        }
    }

    auto getVal = [&](int x, int y, int z) -> double {
        double v = newCCForMarked[x][y][z];

        // Treat uninitialized as 0 to match Mathematica's ConstantArray[0]
        return v;
    };


    // Iterate centers of 3x3x3 neighborhoods: indices 1..(h-2) stepping by 2
    for (int i = 1; i <= hx - 2; i += 2) {
        for (int j = 1; j <= hy - 2; j += 2) {
            for (int k = 1; k <= hz - 2; k += 2) {
                const int mi = (i - 1) / 2;
                const int mj = (j - 1) / 2;
                const int mk = (k - 1) / 2;
                const double tet = getVal(i, j, k);
                // 8 vertices around
                double v1 = getVal(i - 1, j - 1, k - 1);
                double v2 = getVal(i + 1, j - 1, k - 1);
                double v3 = getVal(i - 1, j + 1, k - 1);
                double v4 = getVal(i + 1, j + 1, k - 1);
                double v5 = getVal(i - 1, j - 1, k + 1);
                double v6 = getVal(i + 1, j - 1, k + 1);
                double v7 = getVal(i - 1, j + 1, k + 1);
                double v8 = getVal(i + 1, j + 1, k + 1);
                // 12 edges around
                double e1 = getVal(i - 1, j,     k - 1);
                double e2 = getVal(i + 1, j,     k - 1);
                double e3 = getVal(i,     j - 1, k - 1);
                double e4 = getVal(i,     j + 1, k - 1);
                double e5 = getVal(i - 1, j - 1, k    );
                double e6 = getVal(i + 1, j - 1, k    );
                double e7 = getVal(i - 1, j + 1, k    );
                double e8 = getVal(i + 1, j + 1, k    );
                double e9 = getVal(i - 1, j,     k + 1);
                double e10 = getVal(i + 1, j,    k + 1);
                double e11 = getVal(i,     j - 1, k + 1);
                double e12 = getVal(i,     j + 1, k + 1);
                // 6 faces around
                double f1 = getVal(i - 1, j, k);
                double f2 = getVal(i + 1, j, k);
                double f3 = getVal(i, j - 1, k);
                double f4 = getVal(i, j + 1, k);
                double f5 = getVal(i, j, k - 1);
                double f6 = getVal(i, j, k + 1);
                auto max_all_vertices = std::max({v1,v2,v3,v4,v5,v6,v7,v8});

                if (!sameSign(max_all_vertices,tet)) {
                    marked[mi][mj][mk] = true;
                    continue;
                }
                
                auto max1_2_3_4 = std::max({v1,v2,v3,v4});
                auto max1_5_7_9 = std::max({v1,v3,v5,v7});
                auto max2_6_8_10 = std::max({v2,v4,v6,v8});
                auto max3_5_6_11 = std::max({v1,v2,v5,v6});
                auto max4_7_8_12 = std::max({v3,v4,v7,v8});
                auto max9_10_11_12 = std::max({v5,v6,v7,v8});

                if (!sameSign(max1_5_7_9,f1)) { marked[mi][mj][mk] = true; continue; }
                if (!sameSign(max2_6_8_10,f2)) { marked[mi][mj][mk] = true; continue; }
                if (!sameSign(max3_5_6_11,f3)) { marked[mi][mj][mk] = true; continue; }
                if (!sameSign(max4_7_8_12,f4)) { marked[mi][mj][mk] = true; continue; }
                if (!sameSign(max1_2_3_4,f5))     { marked[mi][mj][mk] = true; continue; }
                if (!sameSign(max9_10_11_12,f6)){ marked[mi][mj][mk] = true; continue; }
                auto max_v1_v3 = std::max(v1, v3);
                auto max_v2_v4 = std::max(v2, v4);
                auto max_v1_v2 = std::max(v1, v2);
                auto max_v3_v4 = std::max(v3, v4);
                auto max_v1_v5 = std::max(v1, v5);
                auto max_v2_v6 = std::max(v2, v6);
                auto max_v3_v7 = std::max(v3, v7);
                auto max_v4_v8 = std::max(v4, v8);
                auto max_v5_v7 = std::max(v5, v7);
                auto max_v6_v8 = std::max(v6, v8);
                auto max_v5_v6 = std::max(v5, v6);
                auto max_v7_v8 = std::max(v7, v8);
                if (sameSign(v1,v3) && !sameSign(max_v1_v3,e1))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v2,v4)&& !sameSign(max_v2_v4,e2))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v1,v2) && !sameSign(max_v1_v2,e3))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v3,v4) && !sameSign(max_v3_v4,e4))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v1,v5) && !sameSign(max_v1_v5,e5))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v2,v6) && !sameSign(max_v2_v6,e6))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v3,v7) && !sameSign(max_v3_v7,e7))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v4,v8) && !sameSign(max_v4_v8,e8))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v5,v7) && !sameSign(max_v5_v7,e9))  { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v6,v8) && !sameSign(max_v6_v8,e10)) { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v5,v6) && !sameSign(max_v5_v6,e11)) { marked[mi][mj][mk] = true; continue; }
                if (sameSign(v7,v8) && !sameSign(max_v7_v8,e12)) { marked[mi][mj][mk] = true; continue; }
            }
        }
    }


    return marked;
}

int runCubicalMode(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] 
                  << " cubical <input_file> <output_file> <output_file2> <adjustment> <dtype> "
                  << "[-a] [--cpp_program <path>] [--topK <value>] "
                  << "[--cavitySkip <n>] [--handleSkip <n>] [--componentSkip <n>]" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = "output/" + std::string(argv[2]);
    std::string output_file2 = "output/" + std::string(argv[3]);
    std::string adjustment_str = argv[4];
    std::string dtype = argv[5];
    
    bool ascii_mode = false;
    std::string cpp_program = "C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\build\\release\\TopoMinCut";
    int topK = 0;
    int cavitySkip = 0;
    int handleSkip = 0;
    int componentSkip = 0;

    // core/neighborhood use floating infinities so callers can check via std::isinf
    double core = std::numeric_limits<double>::infinity();
    double neighborhood = -std::numeric_limits<double>::infinity();
    bool create_3d_array = true;
    
    // Parse optional arguments
    for (int i = 6; i < argc; ++i) {
        if (std::string(argv[i]) == "-a" || std::string(argv[i]) == "--ascii") {
            ascii_mode = true;
        } else if (std::string(argv[i]) == "--cpp_program" && i + 1 < argc) {
            cpp_program = argv[++i];
        } else if (std::string(argv[i]) == "--topK" && i + 1 < argc) {
            topK = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--cavitySkip" && i + 1 < argc) {
            cavitySkip = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--handleSkip" && i + 1 < argc) {
            handleSkip = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "--componentSkip" && i + 1 < argc) {
            componentSkip = std::stoi(argv[++i]);
        }
    }

    double adjustment = std::stod(adjustment_str);

    (void)cpp_program; // TopoMinCut runs in-process; flag kept for backward-compatible CLIs.
    
    core = -(core - adjustment);
    neighborhood = -(neighborhood - adjustment);

    std::cout << "Adjusted Core: " << core << std::endl;
    std::cout << "Adjusted Neighborhood: " << neighborhood << std::endl;
    
    // Read MRC file
    std::cout << "Reading MRC file: " << input_file << std::endl;
    
    std::vector<std::vector<std::vector<double>>> data;
    int dims[3];
    
    try {
        data = MRCReader::readMRCFile(input_file, adjustment);

        dims[0] = data.size();
        dims[1] = data[0].size();
        dims[2] = data[0][0].size();
    } catch (const std::exception& e) {
        std::cerr << "Error reading MRC file: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "Data shape: " << dims[0] << "x" << dims[1] << "x" << dims[2] << std::endl;

    auto old_result = contour3DSmall(data,0.0);
    reportMeshManifoldness(old_result.first, old_result.second, "Old Result");

    int totalVertices = dims[0] * dims[1] * dims[2];

    int xEdgeCount = (dims[0] - 1) * dims[1] * dims[2];
    int yEdgeCount = dims[0] * (dims[1] - 1) * dims[2];
    int zEdgeCount = dims[0] * dims[1] * (dims[2] - 1);
    int totalEdges = xEdgeCount + yEdgeCount + zEdgeCount;

    int xQuadCount = dims[0] * (dims[1] - 1) * (dims[2] - 1);
    int yQuadCount = (dims[0] - 1) * dims[1] * (dims[2] - 1);
    int zQuadCount = (dims[0] - 1) * (dims[1] - 1) * dims[2];
    int totalQuads = xQuadCount + yQuadCount + zQuadCount;

    int totalCubes = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);

    int totalCells = totalVertices + totalEdges + totalQuads + totalCubes;

    std::vector <bool> cellInRange(totalCells, false);

    std::vector <double> cellAlpha(totalCells, 0.0);
    std::vector <double> cellAlphaPlus(totalCells, 0.0);
    std::vector <double> cellAlphaMinus(totalCells, 0.0);

    std::vector <int> cellDimension(totalCells, 0);
    //fill cellDimension with 0, 1, 2, 3 for vertices, edges, quads, and cubes respectively
    for (int i = 0; i < totalVertices; ++i) {
        cellDimension[i] = 0;
    }
    for (int i = totalVertices; i < totalVertices + totalEdges; ++i) {
        cellDimension[i] = 1;
    }
    for (int i = totalVertices + totalEdges; i < totalVertices + totalEdges + totalQuads; ++i) {
        cellDimension[i] = 2;
    }
    for (int i = totalVertices + totalEdges + totalQuads; i < totalCells; ++i) {
        cellDimension[i] = 3;
    }

    //give a vector of vector for parent cell of each cell
    std::vector<std::vector<int>> parentCell(totalCells, std::vector<int>(0));
    //for each cell, reserve 6 spaces for parent cell
    for (int i = 0; i < totalCells; ++i) {
        parentCell[i].reserve(6);
    }

    std::cout << "Identifying inRange and interface cells..." << std::endl;

    // set alpha for vertices
    for (int z = 0; z < dims[2]; ++z) {
        for (int y = 0; y < dims[1]; ++y) {
            for (int x = 0; x < dims[0]; ++x) {
                int id = dims[0] * dims[1] * z + dims[0] * y + x;
                double val = data[x][y][z];
                // Store the value for later comparison
                cellAlpha[id] = val;
            }
        }
    }

    std::cout << "Finished setting up alpha values for verts" << std::endl;

    // Identify inRange edges, quads, and cubes
    for (int z = 0; z < dims[2]; ++z) {
        for (int y = 0; y < dims[1]; ++y) {
            for (int x = 0; x < dims[0]; ++x) {
                int v = dims[0] * dims[1] * z + dims[0] * y + x;
                int vx = v + 1;
                int vy = v + dims[0];
                int vz = v + dims[0] * dims[1];

                int vxy = v + dims[0] + 1;
                int vxz = v + dims[0] * dims[1] + 1;
                int vyz = v + dims[0] * dims[1] + dims[0];
                int vxyz = v + dims[0] * dims[1] + dims[0] + 1;
                
                // Check edges
                if (x < dims[0] - 1) {
                    int edge_id = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices;
                    double edge_val = std::max({cellAlpha[v], cellAlpha[vx]});
                    double edge_alphaMinus = std::min({cellAlpha[v], cellAlpha[vx]});

                    cellAlpha[edge_id] = edge_val;
                    cellAlphaMinus[edge_id] = edge_alphaMinus;

                    parentCell[v].push_back(edge_id);
                    parentCell[vx].push_back(edge_id);
                }
                
                if (y < dims[1] - 1) {
                    int edge_id = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices;
                    double edge_val = std::max({cellAlpha[v], cellAlpha[vy]});
                    double edge_alphaMinus = std::min({cellAlpha[v], cellAlpha[vy]});

                    cellAlpha[edge_id] = edge_val;
                    cellAlphaMinus[edge_id] = edge_alphaMinus;

                    parentCell[v].push_back(edge_id);
                    parentCell[vy].push_back(edge_id);
                }
                if (z < dims[2] - 1) {
                    int edge_id = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                    double edge_val = std::max({cellAlpha[v], cellAlpha[vz]});
                    double edge_alphaMinus = std::min({cellAlpha[v], cellAlpha[vz]});

                    cellAlpha[edge_id] = edge_val;
                    cellAlphaMinus[edge_id] = edge_alphaMinus;

                    parentCell[v].push_back(edge_id);
                    parentCell[vz].push_back(edge_id);
                }
                
                // Check quads
                if (y < dims[1] - 1 && z < dims[2] - 1) {
                    int quad_id = dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices + totalEdges;
                    double quad_val = std::max({cellAlpha[v], cellAlpha[vy], cellAlpha[vz], cellAlpha[vyz]});
                    double quad_alphaMinus = std::min({cellAlpha[v], cellAlpha[vy], cellAlpha[vz], cellAlpha[vyz]});

                    cellAlpha[quad_id] = quad_val;
                    cellAlphaMinus[quad_id] = quad_alphaMinus;

                    int e1 = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x +totalVertices;
                    int e2 = e1 + dims[0] * (dims[1] - 1);
                    int e3 = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                    int e4 = e3 + dims[0];

                    parentCell[e1].push_back(quad_id);
                    parentCell[e2].push_back(quad_id);
                    parentCell[e3].push_back(quad_id);
                    parentCell[e4].push_back(quad_id);
                }
                
                if (x < dims[0] - 1 && z < dims[2] - 1) {
                    int quad_id = xQuadCount + (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    
                    double quad_val = std::max({cellAlpha[v], cellAlpha[vx], cellAlpha[vz], cellAlpha[vxz]});
                    double quad_alphaMinus = std::min({cellAlpha[v], cellAlpha[vx], cellAlpha[vz], cellAlpha[vxz]});

                    cellAlpha[quad_id] = quad_val;
                    cellAlphaMinus[quad_id] = quad_alphaMinus;

                    int e1 = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x +  totalVertices;
                    int e2 = e1 + (dims[0] - 1) * dims[1];
                    int e3 = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                    int e4 = e3 + 1;

                    parentCell[e1].push_back(quad_id);
                    parentCell[e2].push_back(quad_id);
                    parentCell[e3].push_back(quad_id);
                    parentCell[e4].push_back(quad_id);
                }
                if (x < dims[0] - 1 && y < dims[1] - 1) {
                    int quad_id = xQuadCount + yQuadCount + (dims[0] - 1) * (dims[1] - 1) * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;;
                    
                    double quad_val = std::max({cellAlpha[v], cellAlpha[vx], cellAlpha[vy], cellAlpha[vxy]});
                    double quad_alphaMinus = std::min({cellAlpha[v], cellAlpha[vx], cellAlpha[vy], cellAlpha[vxy]});

                    cellAlpha[quad_id] = quad_val;
                    cellAlphaMinus[quad_id] = quad_alphaMinus;

                    int e1 = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices;
                    int e2 = e1 + (dims[0] - 1);
                    int e3 = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices;
                    int e4 = e3 + 1;

                    parentCell[e1].push_back(quad_id);
                    parentCell[e2].push_back(quad_id);
                    parentCell[e3].push_back(quad_id);
                    parentCell[e4].push_back(quad_id);
                }
                
                // Check cubes
                if (x < dims[0] - 1 && y < dims[1] - 1 && z < dims[2] - 1) {
                    int cube_id = z * (dims[0] - 1) * (dims[1] - 1) + y * (dims[0] - 1) + x + totalVertices + totalEdges + totalQuads;

                    double cube_val = std::max({cellAlpha[v], cellAlpha[vx], cellAlpha[vy], cellAlpha[vz], cellAlpha[vxy], cellAlpha[vxz], cellAlpha[vyz], cellAlpha[vxyz]});
                    double cube_alphaMinus = std::min({cellAlpha[v], cellAlpha[vx], cellAlpha[vy], cellAlpha[vz], cellAlpha[vxy], cellAlpha[vxz], cellAlpha[vyz], cellAlpha[vxyz]});
                    
                    cellAlpha[cube_id] = cube_val;
                    cellAlphaMinus[cube_id] = cube_alphaMinus;
                    
                    int q1 = dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices + totalEdges;
                    int q2 = q1 + 1;
                    int q3 = xQuadCount + (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    int q4 = q3 + (dims[0] - 1);
                    int q5 = xQuadCount + yQuadCount + (dims[0] - 1) * (dims[1] - 1) * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    int q6 = q5 + (dims[0] - 1) * (dims[1] - 1);

                    parentCell[q1].push_back(cube_id);
                    parentCell[q2].push_back(cube_id);
                    parentCell[q3].push_back(cube_id);
                    parentCell[q4].push_back(cube_id);
                    parentCell[q5].push_back(cube_id);
                    parentCell[q6].push_back(cube_id);
                    
                }

                
            }
        }
    }

    std::cout << "Finished setting up alphaC values" << std::endl;

    // Process cells by dimension (vertices -> edges -> quads) to enable memoization
    // For each cell, cellAlphaPlus is the max of its parents' cellAlpha (if cube) or cellAlphaPlus (if not cube)

    for (int i = totalCells-1; i >= 0; --i) {
        if (cellDimension[i] == 3) continue;
        
        const auto& parents = parentCell[i];
        if (parents.empty()) cellAlphaPlus[i] = cellAlpha[i];
        
        double maxAlpha = -std::numeric_limits<double>::infinity();
        for (int parentId : parents) {
            double parentAlpha;
            if (cellDimension[parentId] == 3) {
                // Parent is a cube, use its cellAlpha directly
                parentAlpha = cellAlpha[parentId];
            } else {
                // Parent is not a cube, use its cellAlphaPlus (already computed since we process by dimension)
                parentAlpha = cellAlphaPlus[parentId];
            }
            maxAlpha = std::max(maxAlpha, parentAlpha);
        }
        cellAlphaPlus[i] = maxAlpha;
    }
    
    // Now perform interface detection in the correct order
    std::cout << "Performing inRange and interface detection..." << std::endl;

    for (int i = 0; i < totalCells; ++i) {

        if (cellAlpha[i] >= core && cellAlpha[i] <= neighborhood) {
            cellInRange[i] = true;
            continue;
        }

        if (cellDimension[i] == 0) {

            if (cellAlpha[i] < core && cellAlphaPlus[i] >= core) {
                cellInRange[i] = true;
                continue;
            }
        }
        else if (cellDimension[i] == 3) {
            
            if (cellAlphaMinus[i] <= neighborhood && cellAlpha[i] > neighborhood) {
                cellInRange[i] = true;
                continue;
            }
        }

        else{
            if (cellAlpha[i] < core && cellAlphaPlus[i] >= core) {
                cellInRange[i] = true;
                continue;
            }
            if (cellAlphaMinus[i] <= neighborhood && cellAlpha[i] > neighborhood) {
                cellInRange[i] = true;
                continue;
            }
        }

    }

    std::cout << "Finished inRange and interface detection" << std::endl;

    // Now create only the inRange and interface cells
    std::cout << "Creating only inRange and interface cells..." << std::endl;
    //count how many trues are in cellInRange
    int totalIncludedCells = 0;
    for (int i = 0; i < totalCells; ++i) {
        if (cellInRange[i]) {
            totalIncludedCells++;
        }
    }
    std::cout << "Total included cells: " << totalIncludedCells << std::endl;
    
    // Create vertices
    std::vector<Vertex> verts;
    std::vector<int> cellIdToIndex (totalCells, -1);
    
    verts.reserve(totalIncludedCells);

    for (int z = 0; z < dims[2]; ++z) {
        for (int y = 0; y < dims[1]; ++y) {
            for (int x = 0; x < dims[0]; ++x) {
                int id = dims[0] * dims[1] * z + dims[0] * y + x;
                if (cellInRange[id]) {
                    cellIdToIndex[id] = verts.size();
                    verts.emplace_back(verts.size(), cellAlpha[id], 2*x, 2*y, 2*z);
                    
                }
            }
        }
    }
    std::cout << "Finished creating vertices" << std::endl;
    verts.shrink_to_fit();

    // Create edges
    std::vector<Edge*> edges;
    edges.reserve(totalIncludedCells);

    for (int z = 0; z < dims[2]; ++z) {
        for (int y = 0; y < dims[1]; ++y) {
            for (int x = 0; x < dims[0]; ++x) {
                int v = dims[0] * dims[1] * z + dims[0] * y + x;
                int vx = v + 1;
                int vy = v + dims[0];
                int vz = v + dims[0] * dims[1];
                
                if (x < dims[0] - 1) {
                    int edge_id = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices;
                    if (cellInRange[edge_id]) {
                        int v1_idx = cellIdToIndex[v];
                        int v2_idx = cellIdToIndex[vx];
                        if (!(v1_idx == -1 && v2_idx == -1)) {
                            cellIdToIndex[edge_id] = edges.size();

                            //all the non -1 ones should be put into children
                            std::vector<Simplex*> children;
                            if (v1_idx != -1) children.push_back(&verts[v1_idx]);
                            if (v2_idx != -1) children.push_back(&verts[v2_idx]);

                            edges.emplace_back(new Edge(edges.size() + verts.size(), cellAlpha[edge_id], 2*x+1, 2*y, 2*z, children));
                        }
                        else{
                            std::cout << "Edge " << edge_id << " has both vertices as -1 (cell complex error)" << std::endl;
                        }
                    }
                }
                if (y < dims[1] - 1) {
                    int edge_id = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices;
                    if (cellInRange[edge_id]) {
                        int v1_idx = cellIdToIndex[v];
                        int v2_idx = cellIdToIndex[vy];

                        //if they are not both -1
                        if (!(v1_idx == -1 && v2_idx == -1)) {
                            cellIdToIndex[edge_id] = edges.size(); 
                            
                            std::vector<Simplex*> children;
                            if (v1_idx != -1) children.push_back(&verts[v1_idx]);
                            if (v2_idx != -1) children.push_back(&verts[v2_idx]);
                            edges.emplace_back(new Edge(edges.size() + verts.size(), cellAlpha[edge_id], 2*x, 2*y+1, 2*z, children));
                        }
                        else{
                            std::cout << "Edge " << edge_id << " has both vertices as -1 (cell complex error)" << std::endl;
                        }
                    }
                }
                if (z < dims[2] - 1) {
                    int edge_id = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                    if (cellInRange[edge_id]) {
                        int v1_idx = cellIdToIndex[v];
                        int v2_idx = cellIdToIndex[vz];
                        if (!(v1_idx == -1 && v2_idx == -1)) {
                            cellIdToIndex[edge_id] = edges.size();

                            std::vector<Simplex*> children;
                            if (v1_idx != -1) children.push_back(&verts[v1_idx]);
                            if (v2_idx != -1) children.push_back(&verts[v2_idx]);
                            edges.emplace_back(new Edge(edges.size() + verts.size(), cellAlpha[edge_id], 2*x, 2*y, 2*z+1, children));
                        }
                        else{
                            std::cout << "Edge " << edge_id << " has both vertices as -1 (cell complex error)" << std::endl;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Finished creating edges" << std::endl;

    // Create quads
    std::vector<Quad*> quads;
    
    for (int z = 0; z < dims[2]; ++z) {
        for (int y = 0; y < dims[1]; ++y) {
            for (int x = 0; x < dims[0]; ++x) {
                if (y < dims[1] - 1 && z < dims[2] - 1) {
                    int quad_id = dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices + totalEdges;

                    if (cellInRange[quad_id]) {
                        int e1 = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices;
                        int e2 = e1 + dims[0] * (dims[1] - 1);
                        int e3 = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                        int e4 = e3 + dims[0];
                        
                        int e1_idx = cellIdToIndex[e1];
                        int e2_idx = cellIdToIndex[e2];
                        int e3_idx = cellIdToIndex[e3];
                        int e4_idx = cellIdToIndex[e4];

                        if (!(e1_idx == -1 && e2_idx == -1 && e3_idx == -1 && e4_idx == -1)) {
                        
                            cellIdToIndex[quad_id] = quads.size();

                            std::vector<Simplex*> children;
                            if (e1_idx != -1) children.push_back(edges[e1_idx]);
                            if (e2_idx != -1) children.push_back(edges[e2_idx]);
                            if (e3_idx != -1) children.push_back(edges[e3_idx]);
                            if (e4_idx != -1) children.push_back(edges[e4_idx]);

                            quads.emplace_back(new Quad(quads.size() + verts.size() + edges.size(), cellAlpha[quad_id], 2*x, 2*y+1, 2*z+1, children));
                        }
                        else{
                            std::cout << "Face " << quad_id << " has all vertices as -1 (cell complex error)" << std::endl;
                        }

                    }
                }
                if (x < dims[0] - 1 && z < dims[2] - 1) {
                    int quad_id = xQuadCount + (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;

                    if (cellInRange[quad_id]) {
                        int e1 = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices;
                        int e2 = e1 + (dims[0] - 1) * dims[1];
                        int e3 = xEdgeCount + yEdgeCount + dims[0] * dims[1] * z + dims[0] * y + x + totalVertices;
                        int e4 = e3 + 1;
                        
                        int e1_idx = cellIdToIndex[e1];
                        int e2_idx = cellIdToIndex[e2];
                        int e3_idx = cellIdToIndex[e3];
                        int e4_idx = cellIdToIndex[e4];
                        
                        if (!(e1_idx == -1 && e2_idx == -1 && e3_idx == -1 && e4_idx == -1)) {
                        
                            cellIdToIndex[quad_id] = quads.size();

                            std::vector<Simplex*> children;
                            if (e1_idx != -1) children.push_back(edges[e1_idx]);
                            if (e2_idx != -1) children.push_back(edges[e2_idx]);
                            if (e3_idx != -1) children.push_back(edges[e3_idx]);
                            if (e4_idx != -1) children.push_back(edges[e4_idx]);

                            quads.emplace_back(new Quad(quads.size() + verts.size() + edges.size(), cellAlpha[quad_id], 2*x+1, 2*y, 2*z+1, children));
                        }
                        else{
                            std::cout << "Face " << quad_id << " has all vertices as -1 (cell complex error)" << std::endl;
                        }
                    }
                }
                if (x < dims[0] - 1 && y < dims[1] - 1) {
                    int quad_id = xQuadCount + yQuadCount + (dims[0] - 1) * (dims[1] - 1) * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    if (cellInRange[quad_id]) {
                        int e1 = (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices;
                        int e2 = e1 + (dims[0] - 1);
                        int e3 = xEdgeCount + dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices;
                        int e4 = e3 + 1;
                        
                        int e1_idx = cellIdToIndex[e1];
                        int e2_idx = cellIdToIndex[e2];
                        int e3_idx = cellIdToIndex[e3];
                        int e4_idx = cellIdToIndex[e4];
                        
                        if (!(e1_idx == -1 && e2_idx == -1 && e3_idx == -1 && e4_idx == -1)) {
                        
                            cellIdToIndex[quad_id] = quads.size();

                            std::vector<Simplex*> children;
                            if (e1_idx != -1) children.push_back(edges[e1_idx]);
                            if (e2_idx != -1) children.push_back(edges[e2_idx]);
                            if (e3_idx != -1) children.push_back(edges[e3_idx]);
                            if (e4_idx != -1) children.push_back(edges[e4_idx]);

                            quads.emplace_back(new Quad(quads.size() + verts.size() + edges.size(), cellAlpha[quad_id], 2*x+1, 2*y+1, 2*z, children));
                        }
                        else{
                            std::cout << "Face " << quad_id << " has all vertices as -1 (cell complex error)" << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    std::cout << "Finished creating quads" << std::endl;
    
    // Create cubes
    std::vector<Cube> cubes;
    
    for (int z = 0; z < dims[2] - 1; ++z) {
        for (int y = 0; y < dims[1] - 1; ++y) {
            for (int x = 0; x < dims[0] - 1; ++x) {
                int cube_id = z * (dims[0] - 1) * (dims[1] - 1) + y * (dims[0] - 1) + x + totalVertices + totalEdges + totalQuads;
                
                if (cellInRange[cube_id]) {

                    int q1 = dims[0] * (dims[1] - 1) * z + dims[0] * y + x + totalVertices + totalEdges;
                    int q2 = q1 + 1;
                    int q3 = xQuadCount + (dims[0] - 1) * dims[1] * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    int q4 = q3 + (dims[0] - 1);
                    int q5 = xQuadCount + yQuadCount + (dims[0] - 1) * (dims[1] - 1) * z + (dims[0] - 1) * y + x + totalVertices + totalEdges;
                    int q6 = q5 + (dims[0] - 1) * (dims[1] - 1);
                    
                    int q1_idx = cellIdToIndex[q1];
                    int q2_idx = cellIdToIndex[q2];
                    int q3_idx = cellIdToIndex[q3];
                    int q4_idx = cellIdToIndex[q4];
                    int q5_idx = cellIdToIndex[q5];
                    int q6_idx = cellIdToIndex[q6];

                    if (!(q1_idx == -1 && q2_idx == -1 && q3_idx == -1 && q4_idx == -1 && q5_idx == -1 && q6_idx == -1)) {

                        std::vector<Simplex*> children;
                        if (q1_idx != -1) children.push_back(quads[q1_idx]);
                        if (q2_idx != -1) children.push_back(quads[q2_idx]);
                        if (q3_idx != -1) children.push_back(quads[q3_idx]);
                        if (q4_idx != -1) children.push_back(quads[q4_idx]);
                        if (q5_idx != -1) children.push_back(quads[q5_idx]);
                        if (q6_idx != -1) children.push_back(quads[q6_idx]);
                    
                        cubes.emplace_back(cubes.size() + verts.size() + edges.size() + quads.size(), cellAlpha[cube_id],2*x+1, 2*y+1, 2*z+1,children);
                    }
                    else{
                        std::cout << "Cube " << cube_id << " has all vertices as -1 (cell complex error)" << std::endl;
                        std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
                        std::cout << q1 << std::endl;
                        std::cout << q2 << std::endl;
                        std::cout << q3 << std::endl;
                        std::cout << q4 << std::endl;
                        std::cout << q5 << std::endl;
                        std::cout << q6 << std::endl;
                    }

                }
            }
        }
    }
    
    cellAlpha.clear();
    cellAlpha.shrink_to_fit();
    cellAlphaMinus.clear();
    cellAlphaMinus.shrink_to_fit();
    cellAlphaPlus.clear();
    cellAlphaPlus.shrink_to_fit();

    std::cout << "Finished creating cubes" << std::endl;

    // Find quads with only one incident cube (memory efficient)
    std::vector<Quad*> single_cube_quads;
    single_cube_quads.reserve(quads.size() / 6); // Rough estimate
    for (Quad* quad : quads) {
        if (quad == nullptr) continue;
        if (quad->incident_cubes.size() == 1) {
            single_cube_quads.push_back(quad);
        }
    }
    
    // Create a new cube incident to all these quads
    /*
    if (!single_cube_quads.empty()) {
        int new_cube_id = cubes.size();
        
        Quad* first_quad = single_cube_quads[0];

        std::vector<Simplex*> children;

        //0,0,0 is place holder for the centroid of this cube
        Cube new_cube(new_cube_id,std::numeric_limits<double>::infinity(), 0,0,0, children);

        new_cube.val = std::numeric_limits<double>::infinity();
        cubes.push_back(new_cube);
        
        cubes.back().faces.clear();

        // Add faces to the new cube
        for (Quad* quad : single_cube_quads) {
            quad->incident_cubes.push_back(&cubes.back());
            cubes.back().faces.push_back(quad);  // Add to the actual cube in the vector
        }
        
        std::cout << "New cube faces: " << cubes.back().faces.size() << std::endl;
        
        std::cout << "Created new cube incident to " << cubes.back().faces.size() << " quads with single incident cube" << std::endl;
    }
    else{
        std::cout << "No quads with single incident cube found" << std::endl;
    }
    */

    // Create all simplices vector (memory efficient)
    int totalSimplices = verts.size() + edges.size() + quads.size() + cubes.size();
    std::vector<Simplex*> allSimplices;
    allSimplices.reserve(totalSimplices);
    
    for (auto& vertex : verts) {
        allSimplices.push_back(&vertex);
    }
    for (Edge* edge : edges) {
        if (edge != nullptr) {
            allSimplices.push_back(edge);
        }
    }
    for (Quad* quad : quads) {
        if (quad != nullptr) {
            allSimplices.push_back(quad);
        }
    }
    for (auto& cube : cubes) {
        allSimplices.push_back(&cube);
    }
    
    // Optimize memory for all simplices
    std::cout << "Optimizing memory for all simplices..." << std::endl;
    for (auto& vertex : verts) {
        vertex.optimizeMemory();
    }
    for (Edge* edge : edges) {
        if (edge != nullptr) {
            edge->optimizeMemory();
        }
    }
    for (Quad* quad : quads) {
        if (quad != nullptr) {
            quad->optimizeMemory();
        }
    }
    for (auto& cube : cubes) {
        cube.optimizeMemory();
    }
    
    std::cout << "Before Sorting..." << std::endl;
    // Sort all simplices
    std::sort(allSimplices.begin(), allSimplices.end(), 
              [](const Simplex* a, const Simplex* b) { return *a < *b; });
    
    // Create sorted indices (memory efficient)
    std::vector<int> sortedIndices;
    sortedIndices.reserve(totalSimplices);
    sortedIndices.resize(totalSimplices, -1);
    
    std::cout << "Finished Sorting..." << std::endl;
    
    for (size_t i = 0; i < allSimplices.size(); ++i) {
        Simplex* currSimplex = allSimplices[i];
        if(currSimplex->id < 0) {
            std::cout << "ERROR: Simplice " << i << " has id " << currSimplex->id << std::endl;
        }
        if (currSimplex->id >= totalSimplices) {
            std::cout << "ERROR: Simplice " << i << " has id " << currSimplex->id << " which is greater than total simplices " << totalSimplices << std::endl;
        }
        sortedIndices[currSimplex->id] = i;
    }

    std::cout <<"Total size of simplices :   " << allSimplices.size() << std::endl;
    std::cout << "Finished creating sorted indices" << std::endl;

    // Build children-per-column boundary structure from faces relationships (memory efficient)
    std::vector<std::vector<int>> childrenByColumn;
    childrenByColumn.reserve(allSimplices.size());
    childrenByColumn.resize(allSimplices.size());
    
    for (size_t col = 0; col < allSimplices.size(); ++col) {
        const auto& faces = allSimplices[col]->getFaces();
        auto& out = childrenByColumn[col];
        out.reserve(faces.size());
        for (Simplex* f : faces) {
            int sortedID = sortedIndices[f->id];
            if (sortedID >= 0) out.push_back(sortedID);
        }
        
    }

    
    // Free single_cube_quads as it's no longer needed
    single_cube_quads.clear();
    single_cube_quads.shrink_to_fit();

    // Filter simplices by inRange ∪ interface using core/neighborhood
    const double alpha1 = core;
    const double alpha2 = neighborhood;

    //testing things
    FilterResult filt = filterInterfaceAndInRange(allSimplices, childrenByColumn, double(-4.0), double(4.0));

    // Free memory from original data structures that are no longer needed
    std::cout << "Freeing memory from original data structures..." << std::endl;
    
    // Free childrenByColumn as it's been processed
    childrenByColumn.clear();
    childrenByColumn.shrink_to_fit();
    
    // Free sortedIndices as they're no longer needed after filtering
    sortedIndices.clear();
    sortedIndices.shrink_to_fit();
    
    // Free edges and quads vectors as they're now in allSimplices
    edges.clear();
    edges.shrink_to_fit();
    quads.clear();
    quads.shrink_to_fit();

    
    // Persist filtered inputs for TopoMinCut (ASCII mode expected by your pipeline)
    {
        std::ofstream f(output_file);
        f << filt.keptSimplices.size() << std::endl;
        for (size_t col = 0; col < filt.keptSimplices.size(); ++col) {
            f << filt.keptSimplices[col]->dimension;

            if (filt.keptSimplices[col]->dimension!= 0 && filt.keptChildrenByCol[col].size()==0){
                std::cout << "ERROR: Simplice " << col << " has dimension " << filt.keptSimplices[col]->dimension << " and no children" << std::endl;
            }

            for (int row : filt.keptChildrenByCol[col]) {
                f << " " << row;
            }
            f << std::endl;
        }
    }
    {
        std::ofstream f2(output_file2);
        f2 << filt.keptSimplices.size() << std::endl;
        for (auto* s : filt.keptSimplices) {
            f2 << s->val << std::endl;
        }
    }
    
    std::vector<int> newToOld = filt.newToOld;
    // Note: marching cubes on original data removed since data was cleared for memory efficiency
    // If needed, re-read data or use a different approach

    
    std::cout << "Skipping marching cubes on original data for memory efficiency" << std::endl;

    std::cout << "Running TopoMinCut pipeline (in-process, Eigen inputs)..." << std::endl;
    const size_t n_bm = filt.keptSimplices.size();
    std::vector<Eigen::Triplet<int>> bm_triplets;
    bm_triplets.reserve(n_bm * 8);
    for (size_t col = 0; col < n_bm; ++col) {
        for (int row : filt.keptChildrenByCol[col]) {
            bm_triplets.emplace_back(static_cast<int>(row), static_cast<int>(col), 1);
        }
    }
    Eigen::SparseMatrix<int> bm_sparse(static_cast<int>(n_bm), static_cast<int>(n_bm));
    bm_sparse.setFromTriplets(bm_triplets.begin(), bm_triplets.end());

    std::vector<int> bm_dims(n_bm);
    std::vector<double> bm_alphas(n_bm);
    for (size_t c = 0; c < n_bm; ++c) {
        bm_dims[c] = filt.keptSimplices[c]->dimension;
        bm_alphas[c] = filt.keptSimplices[c]->val;
    }

    topomincut::RunParams tm_params;
    tm_params.topK = topK;
    tm_params.core = core;
    tm_params.neighborhood = neighborhood;
    tm_params.cavitySkip = cavitySkip;
    tm_params.handleSkip = handleSkip;
    tm_params.componentSkip = componentSkip;
    topomincut::RunOutputs tm_out;
    const int tm_rc = topomincut::runFromEigenSparse(bm_sparse, bm_dims, bm_alphas, tm_params, &tm_out);
    if (tm_rc != 0) {
        std::cerr << "TopoMinCut pipeline failed with exit code " << tm_rc << std::endl;
        return 1;
    }

    std::vector<double> new_alpha_values = tm_out.alphas_updated;
    /*
    std::vector<int> birth_generating_sets = read_cell_indices("C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\output\\iset_birth_generating_sets.txt");
    std::vector<int> death_generating_sets = read_cell_indices("C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\output\\iset_death_generating_sets.txt");
    */
    //instead of reading a file, I want to read all files in a directory and store them in a vector of vectors
    //so do a for loop
  
    /*
    std::vector<std::vector<int>> generators_all;
    for (const auto& entry : std::filesystem::directory_iterator("C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\output\\generators")) {
        if (entry.is_regular_file()) {
            std::vector<int> generators = read_cell_indices(entry.path().string());
            generators_all.push_back(generators);
        }
    }
    */
    std::vector<int> remaining_negative = tm_out.remaining_negative;
    std::vector<int> protected_indices = tm_out.protected_indices;
    //std::vector<int> new_generators = read_cell_indices("C:\\Users\\l.rong\\Desktop\\Research\\TopoMinCut\\output\\new_generators.txt");


    
    //clear output/embedded_generators/ directory
    //std::filesystem::remove_all("output/embedded_generators/");
    //std::filesystem::create_directory("output/embedded_generators/");
    
    //now instead of a single file do a for loop to write each generator to a separate file in a directory

    // Output in OBJ and PLY format for MeshLab
    {
        std::cout << "Writing remaining negative to OBJ and PLY files..." << std::endl;
        
        // Ensure output directory exists
        std::filesystem::create_directories("output");
        
        // First pass: collect all unique vertices
        std::map<std::tuple<int, int, int>, int> vertexMap;
        std::vector<std::tuple<int, int, int>> vertices;
        std::vector<std::pair<int, int>> edges;
        std::vector<std::vector<int>> faces; // Store face vertex indices
        int vertexIndex = 0; // OBJ uses 1-based indexing, PLY uses 0-based
        
        // Collect vertices from dimension 0 simplices
        for (int id : remaining_negative) {
            if (filt.keptSimplices[id]->dimension == 0) {
                int x = filt.keptSimplices[id]->x / 2;
                int y = filt.keptSimplices[id]->y / 2;
                int z = filt.keptSimplices[id]->z / 2;
                std::tuple<int, int, int> vertex = std::make_tuple(x, y, z);
                
                if (vertexMap.find(vertex) == vertexMap.end()) {
                    vertexMap[vertex] = vertexIndex;
                    vertices.push_back(vertex);
                    vertexIndex++;
                }
            }
        }
        
        // Collect vertices from edges (fixed: check which coordinate is odd for edge direction)
        for (int id : remaining_negative) {
            if (filt.keptSimplices[id]->dimension == 1) {
                int x = filt.keptSimplices[id]->x;
                int y = filt.keptSimplices[id]->y;
                int z = filt.keptSimplices[id]->z;
                
                // For edges: one coordinate is odd (the direction), others are even
                // Edge in x direction: (2*x+1, 2*y, 2*z) -> vertices at (x, y, z) and (x+1, y, z)
                // Edge in y direction: (2*x, 2*y+1, 2*z) -> vertices at (x, y, z) and (x, y+1, z)
                // Edge in z direction: (2*x, 2*y, 2*z+1) -> vertices at (x, y, z) and (x, y, z+1)
                std::tuple<int, int, int> v1, v2;
                if (x % 2 == 1) {
                    // Edge in x direction (x is odd, y and z are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple((x + 1) / 2, y / 2, z / 2);
                } else if (y % 2 == 1) {
                    // Edge in y direction (y is odd, x and z are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple(x / 2, (y + 1) / 2, z / 2);
                } else if (z % 2 == 1) {
                    // Edge in z direction (z is odd, x and y are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple(x / 2, y / 2, (z + 1) / 2);
                } else {
                    std::cout << "ERROR: Edge " << id << " has no odd dimension" << std::endl;
                    continue;
                }
                
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
                //create edge
                edges.push_back(std::make_pair(vertexMap[v1],vertexMap[v2]));
            }
        }
        
        // Process quads and collect any missing vertices, then create faces
        for (int id : remaining_negative) {
            if (filt.keptSimplices[id]->dimension == 2) {
                std::vector<std::tuple<int, int, int>> quadVertices;
                
                // Extract quad vertices based on which dimension is even
                if (filt.keptSimplices[id]->x % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                } else if (filt.keptSimplices[id]->y % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                } else if (filt.keptSimplices[id]->z % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                } else {
                    std::cout << "ERROR: Quad " << id << " has no even dimension" << std::endl;
                    continue;
                }
                
                // Add any missing vertices and get indices
                std::vector<int> indices;
                for (const auto& v : quadVertices) {
                    if (vertexMap.find(v) == vertexMap.end()) {
                        vertexMap[v] = vertexIndex;
                        vertices.push_back(v);
                        vertexIndex++;
                    }
                    indices.push_back(vertexMap[v]);
                }
                
                // Store quad as a single face with 4 indices
                if (indices.size() == 4) {
                    faces.push_back(indices);
                }
            }
        }
        
        // Write PLY file: ASCII format (0-based indexing)
        std::ofstream plyFile("output/output_skeleton_edges.ply");
        plyFile << "ply\n";
        plyFile << "format ascii 1.0\n";
        plyFile << "comment TopoCubical2 output\n";
        plyFile << "element vertex " << vertices.size() << "\n";
        plyFile << "property float x\n";
        plyFile << "property float y\n";
        plyFile << "property float z\n";
        plyFile << "element edge " << edges.size() << "\n";
        plyFile << "property int vertex1\n";
        plyFile << "property int vertex2\n";
        plyFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }
        
        // Write edges (0-based)
        for (const auto& e : edges) {
            plyFile << e.first << " " << e.second << "\n";
        }

        plyFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << faces.size() 
                  << " faces to output/output_skeleton_edges.ply" << std::endl;

        std::ofstream plyFacesFile("output/output_skeleton_faces.ply");
        plyFacesFile << "ply\n";
        plyFacesFile << "format ascii 1.0\n";
        plyFacesFile << "comment TopoCubical2 output\n";
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
                plyFacesFile << " " << (idx); // Convert from 1-based to 0-based
            }
            plyFacesFile << "\n";
        }
        plyFacesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << faces.size() 
                  << " faces to output/output_skeleton_faces.ply" << std::endl;
    }
    
    // Output in OBJ and PLY format for MeshLab
    {
        std::cout << "Writing remaining negative PLY files..." << std::endl;
        
        // Ensure output directory exists
        std::filesystem::create_directories("output");
        
        // First pass: collect all unique vertices
        std::map<std::tuple<int, int, int>, int> vertexMap;
        std::vector<std::tuple<int, int, int>> vertices;
        std::vector<std::pair<int, int>> edges;
        std::vector<std::vector<int>> faces; // Store face vertex indices
        int vertexIndex = 0; //PLY uses 0-based
        
        // Collect vertices from dimension 0 simplices
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 0) {
                int x = filt.keptSimplices[id]->x / 2;
                int y = filt.keptSimplices[id]->y / 2;
                int z = filt.keptSimplices[id]->z / 2;
                std::tuple<int, int, int> vertex = std::make_tuple(x, y, z);
                
                if (vertexMap.find(vertex) == vertexMap.end()) {
                    vertexMap[vertex] = vertexIndex;
                    vertices.push_back(vertex);
                    vertexIndex++;
                }
            }
        }
        
        // Collect vertices from edges (fixed: check which coordinate is odd for edge direction)
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 1) {
                int x = filt.keptSimplices[id]->x;
                int y = filt.keptSimplices[id]->y;
                int z = filt.keptSimplices[id]->z;
                
                // For edges: one coordinate is odd (the direction), others are even
                // Edge in x direction: (2*x+1, 2*y, 2*z) -> vertices at (x, y, z) and (x+1, y, z)
                // Edge in y direction: (2*x, 2*y+1, 2*z) -> vertices at (x, y, z) and (x, y+1, z)
                // Edge in z direction: (2*x, 2*y, 2*z+1) -> vertices at (x, y, z) and (x, y, z+1)
                std::tuple<int, int, int> v1, v2;
                if (x % 2 == 1) {
                    // Edge in x direction (x is odd, y and z are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple((x + 1) / 2, y / 2, z / 2);
                } else if (y % 2 == 1) {
                    // Edge in y direction (y is odd, x and z are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple(x / 2, (y + 1) / 2, z / 2);
                } else if (z % 2 == 1) {
                    // Edge in z direction (z is odd, x and y are even)
                    v1 = std::make_tuple(x / 2, y / 2, z / 2);
                    v2 = std::make_tuple(x / 2, y / 2, (z + 1) / 2);
                } else {
                    std::cout << "ERROR: Edge " << id << " has no odd dimension" << std::endl;
                    continue;
                }
                
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
                //create edge
                edges.push_back(std::make_pair(vertexMap[v1],vertexMap[v2]));
            }
        }
        
        // Process quads and collect any missing vertices, then create faces
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 2) {
                std::vector<std::tuple<int, int, int>> quadVertices;
                
                // Extract quad vertices based on which dimension is even
                if (filt.keptSimplices[id]->x % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                } else if (filt.keptSimplices[id]->y % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        (filt.keptSimplices[id]->z + 1) / 2));
                } else if (filt.keptSimplices[id]->z % 2 == 0) {
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        filt.keptSimplices[id]->y / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        (filt.keptSimplices[id]->x + 1) / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                    quadVertices.push_back(std::make_tuple(
                        filt.keptSimplices[id]->x / 2,
                        (filt.keptSimplices[id]->y + 1) / 2,
                        filt.keptSimplices[id]->z / 2));
                } else {
                    std::cout << "ERROR: Quad " << id << " has no even dimension" << std::endl;
                    continue;
                }
                
                // Add any missing vertices and get indices
                std::vector<int> indices;
                for (const auto& v : quadVertices) {
                    if (vertexMap.find(v) == vertexMap.end()) {
                        vertexMap[v] = vertexIndex;
                        vertices.push_back(v);
                        vertexIndex++;
                    }
                    indices.push_back(vertexMap[v]);
                }
                
                // Store quad as a single face with 4 indices
                if (indices.size() == 4) {
                    faces.push_back(indices);
                }
            }
        }
        
        // Write PLY file: ASCII format (0-based indexing)
        std::ofstream plyLinesFile("output/input_skeleton_edges.ply");
        plyLinesFile << "ply\n";
        plyLinesFile << "format ascii 1.0\n";
        plyLinesFile << "comment TopoCubical2 output\n";
        plyLinesFile << "element vertex " << vertices.size() << "\n";
        plyLinesFile << "property float x\n";
        plyLinesFile << "property float y\n";
        plyLinesFile << "property float z\n";
        plyLinesFile << "element edge " << edges.size() << "\n";
        plyLinesFile << "property int vertex1\n";
        plyLinesFile << "property int vertex2\n";
        plyLinesFile << "end_header\n";
        
        // Write vertices (0-based)
        for (const auto& v : vertices) {
            plyLinesFile << std::get<0>(v) << " " << std::get<1>(v) << " " << std::get<2>(v) << "\n";
        }

        
        // Write edges (0-based)
        for (const auto& e : edges) {
            plyLinesFile << e.first << " " << e.second  << " \n";
        }

        plyLinesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << edges.size() 
                  << " edges to output/input_skeleton_edges.ply" << std::endl;

        
        std::ofstream plyFacesFile("output/input_skeleton_faces.ply");
        plyFacesFile << "ply\n";
        plyFacesFile << "format ascii 1.0\n";
        plyFacesFile << "comment TopoCubical2 output\n";
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
                plyFacesFile << " " << (idx); // Convert from 1-based to 0-based
            }
            plyFacesFile << "\n";
        }
        plyFacesFile.close();
        std::cout << "Wrote " << vertices.size() << " vertices and " << faces.size() 
                  << " faces to output/input_skeleton_faces.ply" << std::endl;
    }
    /*
     //output remaining_negative.txt to mathematica visual
     std::cout << "Writing new generators..." << std::endl;
     // Output inRange cells
     {
         std::ofstream ofs("output/new_generators.wl");
         ofs << "(* Remaining negative: {vertices, edges, quads} *)\n";
         ofs << "{\n";
         // Vertices
         ofs << "  (* Vertices *)\n  {";
         bool first = true;
 
         //loop through all generators_embedded and for all the 0 dimension simplices, write the vertex coordinates
         for (int id : new_generators) {
             if (filt.keptSimplices[id]->dimension == 0) {
                 if (!first) ofs << ", ";
                 ofs << "{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}";
                 first = false;
             }
         }
         ofs << "},\n";
         // Edges (as pairs of vertices)
         ofs << "  (* Edges *)\n  {";
         first = true;
 
         //loop through all generators_embedded and for all the 1 dimension simplices, write the edge coordinates
         for (int id : new_generators) {
             if (filt.keptSimplices[id]->dimension == 1) {
                 if (!first) ofs << ", ";
                 ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
                
                 first = false;
             }
         }
         ofs << "},\n";
         // Quads (as 4 vertices)
         ofs << "  (* Quads *)\n  {";
         first = true;
         //loop through all generators_embedded and for all the 2 dimension simplices, write the quad coordinates
         for (int id : new_generators) {
             if (filt.keptSimplices[id]->dimension == 2) {
                 if (!first) ofs << ", ";
                 //there's going to be 1 even number amongst xyz, that's the one that's going to be the quad direction
                 if (filt.keptSimplices[id]->x % 2 == 0) {
                     ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}, {" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
                 } else if (filt.keptSimplices[id]->y % 2 == 0) {
                     ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}, {" << (filt.keptSimplices[id]->x)/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
                 } else if (filt.keptSimplices[id]->z % 2 == 0) {
                     ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2 << ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}}";
                 } else {
                     std::cout << "ERROR: Quad " << id << " has no even dimension" << std::endl;
                 }
                 first = false;
             }
         }
         ofs << "},\n";
         ofs << "}\n";
         ofs.close();
         std::cout << "Wrote remaining negative to output/new_generators.wl" << std::endl;
     }
    
    //output remaining_negative.txt to mathematica visual
    std::cout << "Writing protected indices to Mathematica files..." << std::endl;
    // Output inRange cells
    {
        std::ofstream ofs("output/protected_indices.wl");
        ofs << "(* protected_indices: {vertices, edges, quads} *)\n";
        ofs << "{\n";
        // Vertices
        ofs << "  (* Vertices *)\n  {";
        bool first = true;

        //loop through all generators_embedded and for all the 0 dimension simplices, write the vertex coordinates
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 0) {
                if (!first) ofs << ", ";
                ofs << "{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}";
                first = false;
            }
        }
        ofs << "},\n";
        // Edges (as pairs of vertices)
        ofs << "  (* Edges *)\n  {";
        first = true;

        //loop through all generators_embedded and for all the 1 dimension simplices, write the edge coordinates
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 1) {
                if (!first) ofs << ", ";
                ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
               
                first = false;
            }
        }
        ofs << "},\n";
        // Quads (as 4 vertices)
        ofs << "  (* Quads *)\n  {";
        first = true;
        //loop through all generators_embedded and for all the 2 dimension simplices, write the quad coordinates
        for (int id : protected_indices) {
            if (filt.keptSimplices[id]->dimension == 2) {
                if (!first) ofs << ", ";
                //there's going to be 1 even number amongst xyz, that's the one that's going to be the quad direction
                if (filt.keptSimplices[id]->x % 2 == 0) {
                    ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}, {" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
                } else if (filt.keptSimplices[id]->y % 2 == 0) {
                    ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}, {" << (filt.keptSimplices[id]->x)/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << (filt.keptSimplices[id]->z+1)/2 << "}}";
                } else if (filt.keptSimplices[id]->z % 2 == 0) {
                    ofs << "{{" << filt.keptSimplices[id]->x/2 << ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << filt.keptSimplices[id]->y/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << (filt.keptSimplices[id]->x+1)/2<< ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}, {" << filt.keptSimplices[id]->x/2 << ", " << (filt.keptSimplices[id]->y+1)/2 << ", " << filt.keptSimplices[id]->z/2 << "}}";
                } else {
                    std::cout << "ERROR: Quad " << id << " has no even dimension" << std::endl;
                }
                first = false;
            }
        }
        ofs << "},\n";
        ofs << "}\n";
        ofs.close();
        std::cout << "Wrote remaining negative to output/protected_indices.wl" << std::endl;
    }
    */
    // Create toHalfVerts array with dimensions = original dims
    std::vector<std::vector<std::vector<double>>> toHalfVerts(dims[0]);
    for (int i = 0; i < dims[0]; ++i) {
        toHalfVerts[i].resize(dims[1]);
        for (int j = 0; j < dims[1]; ++j) {
            toHalfVerts[i][j].resize(dims[2], 0.0);
        }
    }
    
    // Update simplex values with new alpha values
    for (size_t j = 0; j < allSimplices.size() && j < new_alpha_values.size(); ++j) {

        int i = newToOld[j];

        double old_val = allSimplices[i]->val;
        double new_val = new_alpha_values[i];
        if (!sameSign(old_val, new_val)) {
            allSimplices[i]->val = new_val;
            
            if (allSimplices[i]->dimension == 0) {
                //find the vertex actual coords
                Vertex* vertex = static_cast<Vertex*>(allSimplices[i]);
                int x = static_cast<int>(std::round(vertex->x ))/2;
                int y = static_cast<int>(std::round(vertex->y ))/2;
                int z = static_cast<int>(std::round(vertex->z ))/2;
                if (x >= 0 && x < dims[0] && y >= 0 && y < dims[1] && z >= 0 && z < dims[2]) {
                    toHalfVerts[x][y][z] = 1;
                    
                    //also set toHalfVerts to 1 for the vertex in the high-resolution array
                } else {
                    std::cout << "ERROR: Vertex coordinates (" << x << "," << y << "," << z 
                                << ") out of bounds for toHalfVerts array with dimensions [" 
                                << dims[0] << "," << dims[1] << "," << dims[2] << "]" << std::endl;
                }
            }
        }
        
    }

    filt.deleteAll();

    new_alpha_values.clear();
    new_alpha_values.shrink_to_fit();
    newToOld.clear();
    newToOld.shrink_to_fit();
    remaining_negative.clear();
    remaining_negative.shrink_to_fit();

    std::cout << "Finished updating simplex values with new alpha values" << std::endl;
    // Create the high-resolution 3D array first
    auto newCCForMarked = createHighResArray(allSimplices, data, dims);

    // Pad newCCForMarked array with a border of +infinity values on all sides
    
        int nx = static_cast<int>(newCCForMarked.size());
        int ny = static_cast<int>(newCCForMarked[0].size());
        int nz = static_cast<int>(newCCForMarked[0][0].size());

        double maxVal = std::numeric_limits<double>::lowest();
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    maxVal = std::max(maxVal, newCCForMarked[i][j][k]);
                }
            }
        }

        double largeNumber = maxVal + 100;
        int px = nx + 4;
        int py = ny + 4;
        int pz = nz + 4;

        std::vector<std::vector<std::vector<double>>> paddedCC(px);
        for (int i = 0; i < px; ++i) {
            paddedCC[i].resize(py);
            for (int j = 0; j < py; ++j) {
                paddedCC[i][j].resize(pz, largeNumber);
            }
        }
        // Copy data into interior
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    paddedCC[i+2][j+2][k+2] = newCCForMarked[i][j][k];
                }
            }
        }
    
    // Clear newCCForMarked to free memory (no longer needed after padding)
    newCCForMarked.clear();
    newCCForMarked.shrink_to_fit();

    allSimplices.clear();
    allSimplices.shrink_to_fit();

    // Clean up memory
    for (Edge* edge : edges) {
        delete edge;
    }
    for (Quad* quad : quads) {
        delete quad;
    }
    verts.clear();
    edges.clear();
    quads.clear();
    cubes.clear();
    data.clear();

    verts.shrink_to_fit();
    edges.shrink_to_fit();
    quads.shrink_to_fit();
    cubes.shrink_to_fit();
    data.shrink_to_fit();
    std::cout << "Freed newCCForMarked memory" << std::endl;
    
    // Convert high-res array to cell complex structure
    // Use paddedCC if you want padding, otherwise use newCCForMarked

    auto [cellComplex, vertNeedMidPoint] = createCellComplexFromHighRes(paddedCC, toHalfVerts);
    std::cout << "Created cell complex with " << cellComplex.cc[0].size() << " vertices, " 
              << cellComplex.cc[1].size() << " edges, " << cellComplex.cc[2].size() << " quads, "
              << cellComplex.cc[3].size() << " cubes" << std::endl;
    
    // Clear paddedCC to free memory (no longer needed after cell complex creation)
    paddedCC.clear();
    paddedCC.shrink_to_fit();

    std::cout << "Freed paddedCC memory" << std::endl;
    
    int originalNumVertices = cellComplex.cc[0].size();

    auto refinedCellComplex = refineCubesToMixed(cellComplex);
    // Clear cellComplex to free memory (no longer needed after refinement)
    cellComplex.cc.clear();
    cellComplex.vals.clear();
    cellComplex.vertexCoords.clear();
    cellComplex.cc.shrink_to_fit();
    cellComplex.vals.shrink_to_fit();
    cellComplex.vertexCoords.shrink_to_fit();
    std::cout << "Freed cellComplex memory" << std::endl;

    std::cout <<"Refined cell complex has " << refinedCellComplex.vertexCoords.size() << " vertices, " << refinedCellComplex.subs.size() <<" 3 dimensional cells" << std::endl;

    auto finalSurface = contourMixedCells(refinedCellComplex.vertexCoords, refinedCellComplex.vertexVals, refinedCellComplex.subs, originalNumVertices, vertNeedMidPoint);
    // Clear refinedCellComplex to free memory (no longer needed after contouring)
    refinedCellComplex.vertexCoords.clear();
    refinedCellComplex.vertexVals.clear();
    refinedCellComplex.subs.clear();
    refinedCellComplex.vertexCoords.shrink_to_fit();
    refinedCellComplex.vertexVals.shrink_to_fit();
    refinedCellComplex.subs.shrink_to_fit();
    std::cout << "Freed refinedCellComplex memory" << std::endl;
    
    
    // Find special vertices (vertices with val <= 0 that don't have incident edges with val <= 0)
    /*
    std::cout << "Finding special vertices..." << std::endl;
    std::vector<Vertex*> new_special_vertices;
    
    for (Vertex& vertex : verts) {
        if (vertex.val <= 0) {
            bool has_low_val_edge = false;
            
            // Check if any incident edge has new alpha value <= 0
            for (Edge* edge : vertex.incident_edges) {
                if (edge->val <= 0) {
                    has_low_val_edge = true;
                    break;
                }
            }
            
            if (!has_low_val_edge) {
                new_special_vertices.push_back(&vertex);
            }
        }
    }
    
    std::cout << "Found " << new_special_vertices.size() << " special vertices" << std::endl;
    */
    // Create Mathematica visualization files
    std::cout << "Creating Mathematica visualization files..." << std::endl;
    /*
    // 1. Final special elements
    std::ofstream final_special_file("output/final_special_elements.txt");
    if (final_special_file.is_open()) {
        final_special_file << "(* Special elements coordinates after TopoMinCut for Mathematica *)\n";
        
        // Output vertices as points
        final_special_file << "vertices = {\n";
        for (size_t i = 0; i < new_special_vertices.size(); ++i) {
            final_special_file << "  Point[{" << new_special_vertices[i]->x << ", " 
                              << new_special_vertices[i]->y << ", " << new_special_vertices[i]->z << "}]";
            if (i < new_special_vertices.size() - 1) final_special_file << ",";
            final_special_file << "\n";
        }
        final_special_file << "};\n\n";
        
        // Output edges as lines
        final_special_file << "edges = {\n";
        bool has_edges = false;
        for (Simplex* simplex : allSimplices) {
            if (simplex->dimension == 1 && simplex->val <= 0) {
                if (has_edges) final_special_file << ",\n";
                Edge* edge = static_cast<Edge*>(simplex);
                final_special_file << "  Line[{\n";
                final_special_file << "    {" << edge->p1[0] << ", " << edge->p1[1] << ", " << edge->p1[2] << "},\n";
                final_special_file << "    {" << edge->p2[0] << ", " << edge->p2[1] << ", " << edge->p2[2] << "}\n";
                final_special_file << "  }]";
                has_edges = true;
            }
        }
        final_special_file << "\n};\n\n";
        
        // Output faces as polygons
        final_special_file << "faces = {\n";
        bool has_faces = false;
        for (Simplex* simplex : allSimplices) {
            if (simplex->dimension == 2 && simplex->val <= 0) {
                if (has_faces) final_special_file << ",\n";
                Quad* quad = static_cast<Quad*>(simplex);
                final_special_file << "  Polygon[{\n";
                final_special_file << "    {" << quad->p1[0] << ", " << quad->p1[1] << ", " << quad->p1[2] << "},\n";
                final_special_file << "    {" << quad->p2[0] << ", " << quad->p2[1] << ", " << quad->p2[2] << "},\n";
                final_special_file << "    {" << quad->p3[0] << ", " << quad->p3[1] << ", " << quad->p3[2] << "},\n";
                final_special_file << "    {" << quad->p4[0] << ", " << quad->p4[1] << ", " << quad->p4[2] << "}\n";
                final_special_file << "  }]";
                has_faces = true;
            }
        }
        final_special_file << "\n};\n\n";
        
        final_special_file << "(* To visualize all elements, use: *)\n";
        final_special_file << "Graphics3D[{\n";
        final_special_file << "  {Black, PointSize[0.005], vertices},\n";
        final_special_file << "  {Black,Opacity[0.7],Thickness[0.001], edges},\n";
        final_special_file << "  {Black,FaceForm[Opacity[0.1]], EdgeForm[None], faces}\n";
        final_special_file << "}, PlotRange -> All]\n";
        
        final_special_file.close();
        std::cout << "Created final_special_elements.txt" << std::endl;
    }
    
    // 2. Fore generating sets
    std::ofstream fore_file("output/fore_generating_sets.txt");
    if (fore_file.is_open()) {
        fore_file << "(* Fore generating sets visualization for Mathematica *)\n\n";
        
        fore_file << "foreEdges = {\n";
        bool has_fore_edges = false;
        for (int index : generators) {
            if (index < allSimplices.size()) {
                Simplex* simplex = allSimplices[index];
                if (simplex->dimension == 1) { // Include edges
                    if (has_fore_edges) fore_file << ",\n";
                    Edge* edge = static_cast<Edge*>(simplex);
                    fore_file << "  Line[{\n";
                    fore_file << "    {" << edge->p1[0] << ", " << edge->p1[1] << ", " << edge->p1[2] << "},\n";
                    fore_file << "    {" << edge->p2[0] << ", " << edge->p2[1] << ", " << edge->p2[2] << "}\n";
                    fore_file << "  }]";
                    has_fore_edges = true;
                }
            }
        }
        fore_file << "\n};\n\n";

        fore_file << "foreFaces = {\n";
        bool has_fore_faces = false;
        for (int index : generators) {
            if (index < allSimplices.size()) {
                Simplex* simplex = allSimplices[index];
                if (simplex->dimension == 2) { // Include faces
                    if (has_fore_faces) fore_file << ",\n";
                    Quad* quad = static_cast<Quad*>(simplex);
                    fore_file << "  Polygon[{\n";
                    fore_file << "    {" << quad->p1[0] << ", " << quad->p1[1] << ", " << quad->p1[2] << "},\n";
                    fore_file << "    {" << quad->p2[0] << ", " << quad->p2[1] << ", " << quad->p2[2] << "},\n";
                    fore_file << "    {" << quad->p3[0] << ", " << quad->p3[1] << ", " << quad->p3[2] << "},\n";
                    fore_file << "    {" << quad->p4[0] << ", " << quad->p4[1] << ", " << quad->p4[2] << "}\n";
                    fore_file << "  }]";
                    has_fore_faces = true;
                }
            }
        }
        fore_file << "\n};\n\n";
        
        fore_file << "(* To visualize fore generating sets, use: *)\n";
        fore_file << "Graphics3D[{\n";
        fore_file << "  {Yellow, Opacity[0.7], Thickness[0.001], foreEdges},\n";
        fore_file << "  {Yellow, FaceForm[Opacity[0.1]], EdgeForm[None], foreFaces}\n";
        fore_file << "}, PlotRange -> All]\n";
        
        fore_file.close();
        std::cout << "Created fore_generating_sets.txt" << std::endl;
    }
    
    // 3. Birth and death cells (if available)
    if (!birth_generating_sets.empty() || !death_generating_sets.empty()) {
        std::ofstream birth_death_file("output/iset_birth_death_cells.txt");
        if (birth_death_file.is_open()) {
            birth_death_file << "(* Birth and death cells visualization for Mathematica *)\n\n";
            
            // Birth points (vertices)
            birth_death_file << "birthPoints = {\n";
            bool has_birth_points = false;
            for (int index : birth_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 0) {
                        if (has_birth_points) birth_death_file << ",\n";
                        Vertex* v = static_cast<Vertex*>(simplex);
                        birth_death_file << "  Point[{" << v->x << ", " << v->y << ", " << v->z << "}]";
                        has_birth_points = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";

            birth_death_file << "birthEdges = {\n";
            bool has_birth_edges = false;
            for (int index : birth_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 1) {
                        if (has_birth_edges) birth_death_file << ",\n";
                        Edge* edge = static_cast<Edge*>(simplex);
                        birth_death_file << "  Line[{\n";
                        birth_death_file << "    {" << edge->p1[0] << ", " << edge->p1[1] << ", " << edge->p1[2] << "},\n";
                        birth_death_file << "    {" << edge->p2[0] << ", " << edge->p2[1] << ", " << edge->p2[2] << "}\n";
                        birth_death_file << "  }]";
                        has_birth_edges = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";
            
            // Birth faces (quads)
            birth_death_file << "birthFaces = {\n";
            bool has_birth_faces = false;
            for (int index : birth_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 2) {
                        if (has_birth_faces) birth_death_file << ",\n";
                        Quad* quad = static_cast<Quad*>(simplex);
                        birth_death_file << "  Polygon[{\n";
                        birth_death_file << "    {" << quad->p1[0] << ", " << quad->p1[1] << ", " << quad->p1[2] << "},\n";
                        birth_death_file << "    {" << quad->p2[0] << ", " << quad->p2[1] << ", " << quad->p2[2] << "},\n";
                        birth_death_file << "    {" << quad->p3[0] << ", " << quad->p3[1] << ", " << quad->p3[2] << "},\n";
                        birth_death_file << "    {" << quad->p4[0] << ", " << quad->p4[1] << ", " << quad->p4[2] << "}\n";
                        birth_death_file << "  }]";
                        has_birth_faces = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";
            
            // Death points (vertices)
            birth_death_file << "deathPoints = {\n";
            bool has_death_points = false;
            for (int index : death_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 0) {
                        if (has_death_points) birth_death_file << ",\n";
                        Vertex* v = static_cast<Vertex*>(simplex);
                        birth_death_file << "  Point[{" << v->x << ", " << v->y << ", " << v->z << "}]";
                        has_death_points = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";

            birth_death_file << "deathEdges = {\n";
            bool has_death_edges = false;
            for (int index : death_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 1) {
                        if (has_death_edges) birth_death_file << ",\n";
                        Edge* edge = static_cast<Edge*>(simplex);
                        birth_death_file << "  Line[{\n";
                        birth_death_file << "    {" << edge->p1[0] << ", " << edge->p1[1] << ", " << edge->p1[2] << "},\n";
                        birth_death_file << "    {" << edge->p2[0] << ", " << edge->p2[1] << ", " << edge->p2[2] << "}\n";
                        birth_death_file << "  }]";
                        has_death_edges = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";
            
            // Death faces (quads)
            birth_death_file << "deathFaces = {\n";
            bool has_death_faces = false;
            for (int index : death_generating_sets) {
                if (index < allSimplices.size()) {
                    Simplex* simplex = allSimplices[index];
                    if (simplex->dimension == 2) {
                        if (has_death_faces) birth_death_file << ",\n";
                        Quad* quad = static_cast<Quad*>(simplex);
                        birth_death_file << "  Polygon[{\n";
                        birth_death_file << "    {" << quad->p1[0] << ", " << quad->p1[1] << ", " << quad->p1[2] << "},\n";
                        birth_death_file << "    {" << quad->p2[0] << ", " << quad->p2[1] << ", " << quad->p2[2] << "},\n";
                        birth_death_file << "    {" << quad->p3[0] << ", " << quad->p3[1] << ", " << quad->p3[2] << "},\n";
                        birth_death_file << "    {" << quad->p4[0] << ", " << quad->p4[1] << ", " << quad->p4[2] << "}\n";
                        birth_death_file << "  }]";
                        has_death_faces = true;
                    }
                }
            }
            birth_death_file << "\n};\n\n";
            
            birth_death_file << "(* To visualize birth and death cells, use: *)\n";
            birth_death_file << "Graphics3D[{\n";
            birth_death_file << "  {Red, PointSize[0.002], birthPoints},\n";
            birth_death_file << "  {Red, Opacity[0.7], Thickness[0.001], birthEdges},\n";
            birth_death_file << "  {Red, FaceForm[Opacity[0.1]], EdgeForm[None], birthFaces},\n";
            birth_death_file << "  {Blue, PointSize[0.002], deathPoints},\n";
            birth_death_file << "  {Blue, Opacity[0.7], Thickness[0.001], deathEdges},\n";
            birth_death_file << "  {Blue, FaceForm[Opacity[0.1]], EdgeForm[None], deathFaces}\n";
            birth_death_file << "}, PlotRange -> All]\n";
            
            birth_death_file.close();
            std::cout << "Created birth_death_cells.txt" << std::endl;
        }
    }
    */


    std::cout << "Processing completed successfully" << std::endl;
    return 0;

}

} // namespace cubical
} // namespace prescribed_topo
