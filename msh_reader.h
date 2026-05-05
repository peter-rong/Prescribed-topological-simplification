#ifndef MSH_READER_H
#define MSH_READER_H

#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <unordered_set>

struct Point3D {
    double x, y, z;
    Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

struct TetElement {
    int v1, v2, v3, v4;  // vertex indices (1-based in Gmsh, we'll convert to 0-based)
    TetElement(int v1 = 0, int v2 = 0, int v3 = 0, int v4 = 0) : v1(v1), v2(v2), v3(v3), v4(v4) {}
};

class MSHReader {
public:
    static std::pair<std::vector<Point3D>, std::vector<TetElement>> readMSHFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open MSH file: " + filename);
        }

        std::vector<Point3D> vertices;
        std::vector<TetElement> tetrahedra;
        
        std::string line;
        bool inNodes = false;
        bool inElements = false;
        int numNodes = 0;
        int nodesRead = 0;
        int numElements = 0;
        int elementsRead = 0;

        while (std::getline(file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            if (line.empty()) continue;

            // Check for section headers
            if (line == "$Nodes" || line == "$NOD") {
                inNodes = true;
                inElements = false;
                std::getline(file, line);
                std::istringstream iss(line);
                iss >> numNodes;
                vertices.reserve(numNodes);
                nodesRead = 0;
                continue;
            }
            
            if (line == "$EndNodes" || line == "$ENDNOD") {
                inNodes = false;
                continue;
            }

            if (line == "$Elements" || line == "$ELM") {
                inElements = true;
                inNodes = false;
                std::getline(file, line);
                std::istringstream iss(line);
                iss >> numElements;
                tetrahedra.reserve(numElements / 4); // Rough estimate, most elements might be tets
                elementsRead = 0;
                continue;
            }

            if (line == "$EndElements" || line == "$ENDELM") {
                inElements = false;
                continue;
            }

            // Read nodes
            if (inNodes && nodesRead < numNodes) {
                std::istringstream iss(line);
                int nodeId;
                double x, y, z;
                iss >> nodeId >> x >> y >> z;
                vertices.emplace_back(x, y, z);
                nodesRead++;
            }

            // Read elements (we only care about tetrahedra, type 4)
            if (inElements && elementsRead < numElements) {
                std::istringstream iss(line);
                int elementId, elementType, numTags;
                iss >> elementId >> elementType >> numTags;
                
                // Skip tags
                for (int i = 0; i < numTags; ++i) {
                    int tag;
                    iss >> tag;
                }

                // Type 4 is tetrahedron in Gmsh format
                if (elementType == 4) {
                    int v1, v2, v3, v4;
                    iss >> v1 >> v2 >> v3 >> v4;
                    // Convert from 1-based to 0-based indexing
                    tetrahedra.emplace_back(v1 - 1, v2 - 1, v3 - 1, v4 - 1);
                }
                
                elementsRead++;
            }
        }

        file.close();
        
        std::cout << "Read " << vertices.size() << " vertices and " << tetrahedra.size() << " tetrahedra" << std::endl;
        
        return {vertices, tetrahedra};
    }

    // Read vertex alpha values: each line contains a double value (alpha value)
    // Vertex index is the line number (0-based)

    static std::vector<double> readVertexAlphas(const std::string& filename, double adjustment) {
        std::vector<double> vertexAlphas;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open vertex alpha values file: " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            try {
                double alpha = std::stod(line) - adjustment;
                vertexAlphas.push_back(alpha);
            } catch (const std::exception& e) {
                // Skip invalid lines
                continue;
            }
        }
        
        file.close();
        std::cout << "Read " << vertexAlphas.size() << " vertex alpha values" << std::endl;
        
        return vertexAlphas;
    }

    // Read sign file: each line contains an integer value (-1, 0, or 1)
    // Sign index is the line number (0-based)
    static std::vector<int> readSignFile(const std::string& filename) {
        std::vector<int> signs;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open sign file: " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            if (line.empty()) continue;
            
            try {
                int sign = std::stoi(line);
                // Validate that sign is -1, 0, or 1
                if (sign == -1 || sign == 0 || sign == 1) {
                    signs.push_back(sign);
                } else {
                    std::cerr << "Warning: Invalid sign value " << sign << " (expected -1, 0, or 1), skipping" << std::endl;
                }
            } catch (const std::exception& e) {
                // Skip invalid lines
                std::cerr << "Warning: Could not parse line as integer: " << line << std::endl;
                continue;
            }
        }
        
        file.close();
        
        return signs;
    }

    // Read tet labels: each line contains an integer label
    // Tet index is the line number (0-based)
    static std::vector<int> readTetLabels(const std::string& filename) {
        std::vector<int> labels;
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open tet labels file: " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            if (line.empty()) continue;
            
            try {
                int label = std::stoi(line);
                labels.push_back(label);
            } catch (const std::exception& e) {
                // Skip invalid lines
                std::cerr << "Warning: Could not parse line as integer: " << line << std::endl;
                continue;
            }
        }
        
        file.close();
        std::cout << "Read " << labels.size() << " tet labels" << std::endl;

        // INSERT_YOUR_CODE
        int plusOneCount = 0;
        int minusOneCount = 0;
        for (int label : labels) {
            if (label == 1)
                ++plusOneCount;
            else if (label == -1)
                ++minusOneCount;
        }
        std::cout << "Count of +1 labels: " << plusOneCount << std::endl;
        std::cout << "Count of -1 labels: " << minusOneCount << std::endl;
        
        return labels;
    }
};

#endif // MSH_READER_H

