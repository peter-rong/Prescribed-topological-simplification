#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

// Data structures
struct Point3D {
    double x, y, z;
    Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

struct Triangle {
    int v1, v2, v3;
    Triangle(int v1 = 0, int v2 = 0, int v3 = 0) : v1(v1), v2(v2), v3(v3) {}
};

struct MCEdge {
    int p1, p2;
    MCEdge(int p1 = 0, int p2 = 0) : p1(p1), p2(p2) {}
    
    bool operator==(const MCEdge& other) const {
        return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
    }
};

// Hash function for MCEdge
struct MCEdgeHash {
    std::size_t operator()(const MCEdge& edge) const {
        int min_p = std::min(edge.p1, edge.p2);
        int max_p = std::max(edge.p1, edge.p2);
        return std::hash<int>()(min_p) ^ (std::hash<int>()(max_p) << 1);
    }
};

// Marching cubes lookup tables
extern std::vector<std::vector<std::vector<std::pair<int, int>>>> cycles;


// Main function declarations
std::pair<std::vector<Point3D>, std::vector<Triangle>> contour3DLarge(
    const std::vector<std::vector<std::vector<double>>>& img,
    double threshold,
    const std::vector<std::vector<std::vector<bool>>>& marked,
    const std::vector<std::vector<std::vector<double>>>& toHalfVerts);

std::pair<std::vector<Point3D>, std::vector<Triangle>> contour3DSmall(
    const std::vector<std::vector<std::vector<double>>>& img,
    double threshold);
// Helper functions
// Forward declarations for hash structures
struct EdgeCoordinates;
struct EdgeCoordinatesHash;

int getHashDiv(int x1, int y1, int z1, int x2, int y2, int z2,
               const std::vector<std::vector<std::vector<double>>>& zeroImg,
               std::vector<Point3D>& vertices,
               std::unordered_map<EdgeCoordinates, int, EdgeCoordinatesHash>& hash);

int getHash(int x1, int y1, int z1, int x2, int y2, int z2,
            const std::vector<std::vector<std::vector<double>>>& zeroImg,
            const std::vector<std::vector<std::vector<double>>>& toHalfVerts,
            std::vector<Point3D>& vertices,
            std::unordered_map<EdgeCoordinates, int, EdgeCoordinatesHash>& hash);

int getHashSmall(int x1, int y1, int z1, int x2, int y2, int z2,
            const std::vector<std::vector<std::vector<double>>>& zeroImg,
            std::vector<Point3D>& vertices,
            std::unordered_map<EdgeCoordinates, int, EdgeCoordinatesHash>& hash);

Point3D calculateCentroid(const std::vector<Point3D>& vertices, const std::vector<int>& cycle);

bool hasAdjacentMarked(int i, int j, int k, 
                      const std::vector<std::vector<std::vector<bool>>>& marked,
                      int markedDimX, int markedDimY, int markedDimZ);

// Compute marked cubes from a high-resolution scalar field (2x-1 grid)
std::vector<std::vector<std::vector<bool>>> computeMarkedFromHighRes2(
    const std::vector<std::vector<std::vector<double>>>& newCCForMarked);

// Utility: flatten 3D index (x,y,z) for a flat vector with dims (nx,ny,nz)
template <typename T>
inline size_t idx3D(int x, int y, int z, int nx, int ny, int nz) {
    return static_cast<size_t>((x * ny + y) * nz + z);
}

// Write a cubic subarray centered at (cx,cy,cz) with given radius to a Mathematica-readable file
// Data is provided as a flat vector of size nx*ny*nz in x-major order
template <typename T>
inline bool writeMathematica3D(const std::vector<std::vector<std::vector<T>>>& data,
                        int nx, int ny, int nz,
                        int cx, int cy, int cz,
                        int radius,
                        const std::string& outPath)
{
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    int x0 = std::max(0, cx - radius);
    int x1 = std::min(nx - 1, cx + radius);
    int y0 = std::max(0, cy - radius);
    int y1 = std::min(ny - 1, cy + radius);
    int z0 = std::max(0, cz - radius);
    int z1 = std::min(nz - 1, cz + radius);

    std::ofstream ofs(outPath);
    if (!ofs) return false;

    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);

    ofs << "{";
    for (int xi = x0; xi <= x1; ++xi) {
        if (xi != x0) ofs << ",";
        ofs << "{";
        for (int yi = y0; yi <= y1; ++yi) {
            if (yi != y0) ofs << ",";
            ofs << "{";
            for (int zi = z0; zi <= z1; ++zi) {
                if (zi != z0) ofs << ",";
                
                ofs << data[xi][yi][zi];
            }
            ofs << "}";
        }
        ofs << "}";
    }
    ofs << "}";
    ofs.close();
    return true;
}

// Overload 2: write from a flat vector in x-major order
template <typename T>
inline bool writeMathematica3D(const std::vector<T>& data,
                        int nx, int ny, int nz,
                        int cx, int cy, int cz,
                        int radius,
                        const std::string& outPath)
{
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    int x0 = std::max(0, cx - radius);
    int x1 = std::min(nx - 1, cx + radius);
    int y0 = std::max(0, cy - radius);
    int y1 = std::min(ny - 1, cy + radius);
    int z0 = std::max(0, cz - radius);
    int z1 = std::min(nz - 1, cz + radius);

    std::ofstream ofs(outPath);
    if (!ofs) return false;

    ofs << std::setprecision(std::numeric_limits<double>::max_digits10);

    ofs << "{";
    for (int xi = x0; xi <= x1; ++xi) {
        if (xi != x0) ofs << ",";
        ofs << "{";
        for (int yi = y0; yi <= y1; ++yi) {
            if (yi != y0) ofs << ",";
            ofs << "{";
            for (int zi = z0; zi <= z1; ++zi) {
                if (zi != z0) ofs << ",";
                ofs << data[idx3D<T>(xi, yi, zi, nx, ny, nz)];
            }
            ofs << "}";
        }
        ofs << "}";
    }
    ofs << "}";
    ofs.close();
    return true;
}


#endif // MARCHING_CUBES_H
