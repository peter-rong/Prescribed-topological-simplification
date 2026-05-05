#pragma once

#include <Eigen/Sparse>
#include <string>
#include <array>

class BoundaryMatrix {
public:
    BoundaryMatrix();
    
    // Load boundary matrix from file
    bool loadFromFile(const std::string& filename);
    
    // Load alpha values from file
    bool loadAlphas(const std::string& filename);

    /// Square sparse boundary matrix (rows/cols = simplex count) plus simplex dimension per column index.
    bool loadFromEigenCopy(Eigen::SparseMatrix<int> sparse_matrix, std::vector<int> dimensions);

    /// Alpha per column; length must match matrix dimension after loadFromEigenCopy / loadFromFile.
    bool loadAlphasVector(const std::vector<double>& alphas);
    
    // Get the sparse matrix representation
    const Eigen::SparseMatrix<int>& getMatrix() const;
    
    // Get dimensions
    int rows() const;
    int cols() const;

    // Get number of non-zeros in a specific row/column
    int getRowNonZeros(int row) const;
    int getColNonZeros(int col) const;

    // Get all nonzero indices in a specific column
    void getColIndices(int col, std::vector<int>& out) const;
    // Get all nonzero indices in a specific row
    void getRowIndices(int row, std::vector<int>& out) const;

    // Get the alpha values
    const std::vector<double>& getAlphas() const;
    int getAlphaSize() const;

    const std::vector<int>& getDimensions() const;

private:
    Eigen::SparseMatrix<int> matrix_;
    std::vector<double> alphas_; // Dynamic array for alpha values
    int matrix_size_; // Size of the matrix (n x n)
    int alpha_size_; // Size of alpha array
    std::vector<int> dimensions_;
}; 