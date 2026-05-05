#include "boundary_matrix.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <cctype>
#include <algorithm>

// Helper function to convert string to lowercase
std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

BoundaryMatrix::BoundaryMatrix() : matrix_size_(0), alpha_size_(0) {}

bool BoundaryMatrix::loadFromEigenCopy(Eigen::SparseMatrix<int> sparse_matrix, std::vector<int> dimensions) {
    if (sparse_matrix.rows() != sparse_matrix.cols()) {
        std::cerr << "Boundary matrix must be square\n";
        return false;
    }
    const int n = static_cast<int>(sparse_matrix.cols());
    if (static_cast<int>(dimensions.size()) != n) {
        std::cerr << "Dimension count must match matrix size\n";
        return false;
    }
    matrix_size_ = n;
    matrix_ = std::move(sparse_matrix);
    dimensions_ = std::move(dimensions);
    alphas_.clear();
    alpha_size_ = 0;
    return true;
}

bool BoundaryMatrix::loadAlphasVector(const std::vector<double>& alphas) {
    if (matrix_size_ <= 0) {
        std::cerr << "loadAlphasVector: matrix not initialized\n";
        return false;
    }
    if (alphas.size() != static_cast<size_t>(matrix_size_)) {
        std::cerr << "Alpha count " << alphas.size() << " != matrix size " << matrix_size_ << "\n";
        return false;
    }
    alphas_ = alphas;
    alpha_size_ = matrix_size_;
    return true;
}

bool BoundaryMatrix::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open matrix file: " << filename << std::endl;
        return false;
    }

    // Read matrix size from first line
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read matrix size" << std::endl;
        return false;
    }

    // Trim line endings like PHAT does
    line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

    try {
        matrix_size_ = std::stoi(line);
        if (matrix_size_ <= 0) {
            std::cerr << "Invalid matrix size: " << matrix_size_ << std::endl;
            return false;
        }
    } catch (const std::exception& e) {
        std::cerr << "Failed to parse matrix size: " << e.what() << std::endl;
        return false;
    }

    // Initialize sparse matrix
    matrix_ = Eigen::SparseMatrix<int>(matrix_size_, matrix_size_);
    std::vector<Eigen::Triplet<int>> triplets;
    dimensions_.resize(matrix_size_);

    // Read each column
    int col = 0;
    while (std::getline(file, line)) {
        // Trim line endings like PHAT does
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        
        // Skip empty lines and comments like PHAT does
        if (line == "" || line[0] == '#') {
            continue;
        }

        if (col >= matrix_size_) {
            std::cerr << "Found more columns than matrix size" << std::endl;
            return false;
        }

        std::istringstream iss(line);
        int dim;

        if (!(iss >> dim)) {
            std::cerr << "Failed to read dimension at column " << col << std::endl;
            return false;
        }

        dimensions_[col] = dim;

        // Read all row indices into a vector first
        std::vector<int> rows;
        int row;
        while (iss.good()) {
            iss >> row;
            if (row < 0 || row >= matrix_size_) {
                std::cerr << "Invalid row index " << row << " at column " << col << std::endl;
                return false;
            }
            rows.push_back(row);
        }

        // Sort the rows like PHAT does
        std::sort(rows.begin(), rows.end());
        
        // Add to triplets
        for (int row : rows) {
            triplets.emplace_back(row, col, 1);
        }
        col++;
    }

    if (col != matrix_size_) {
        std::cerr << "Expected " << matrix_size_ << " columns but found " << col << std::endl;
        return false;
    }

    matrix_.setFromTriplets(triplets.begin(), triplets.end());
    file.close();
    return true;
}

bool BoundaryMatrix::loadAlphas(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open alpha file: " << filename << std::endl;
        return false;
    }

    // Read size from first line
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read alpha size" << std::endl;
        return false;
    }

    try {
        alpha_size_ = std::stoi(line);
        if (alpha_size_ <= 0) {
            std::cerr << "Invalid alpha size: " << alpha_size_ << std::endl;
            return false;
        }
    } catch (const std::exception& e) {
        std::cerr << "Failed to parse alpha size: " << e.what() << std::endl;
        return false;
    }

    // Resize vector to hold all alpha values
    alphas_.resize(alpha_size_);
    // Read alpha values
    int i = 0;
    while (std::getline(file, line) && i < alpha_size_) {
        try {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            
            // Check for infinity
            if (toLower(line) == "infinity" || toLower(line) == "inf") {
                alphas_[i++] = std::numeric_limits<double>::infinity();
            } else {
                double value = std::stod(line);
                if (!std::isfinite(value) && value > 0) { // Allow positive infinity
                    alphas_[i++] = std::numeric_limits<double>::infinity();
                } else if (!std::isfinite(value)) {
                    std::cerr << "Invalid alpha value at index " << i << ": " << value << std::endl;
                    return false;
                } else {
                    alphas_[i++] = value;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Failed to parse alpha value at index " << i << ": " << e.what() << std::endl;
            return false;
        }
    }

    if (i != alpha_size_) {
        std::cerr << "Expected " << alpha_size_ << " alpha values but found " << i << std::endl;
        return false;
    }
    file.close();
    return true;
}

const Eigen::SparseMatrix<int>& BoundaryMatrix::getMatrix() const {
    return matrix_;
}

int BoundaryMatrix::rows() const {
    return matrix_.rows();
}

int BoundaryMatrix::cols() const {
    return matrix_.cols();
}

int BoundaryMatrix::getRowNonZeros(int row) const {
    if (row < 0 || row >= matrix_.rows()) {
        return 0;
    }
    // Count non-zeros in the row by iterating through the inner vector
    int count = 0;
    for (int j = 0; j < matrix_.cols(); ++j) {
        if (matrix_.coeff(row, j) != 0) {
            count++;
        }
    }
    return count;
}

int BoundaryMatrix::getColNonZeros(int col) const {
    if (col < 0 || col >= matrix_.cols()) {
        return 0;
    }
    // Count non-zeros in the column by iterating through the outer vector
    int count = 0;
    for (int i = 0; i < matrix_.rows(); ++i) {
        if (matrix_.coeff(i, col) != 0) {
            count++;
        }
    }
    return count;
}

const std::vector<double>& BoundaryMatrix::getAlphas() const {
    return alphas_;
}

int BoundaryMatrix::getAlphaSize() const {
    return alpha_size_;
}

void BoundaryMatrix::getColIndices(int col, std::vector<int>& out) const {
    out.clear();
    
    if (col < 0 || col >= matrix_.cols()) {
        return;
    }

    // Count non-zeros in this column
    int nnz = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator it(matrix_, col); it; ++it) {
        nnz++;
    }
    
    // Pre-allocate exact size
    out.resize(nnz);
    
    // Fill directly using iterator
    int idx = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator it(matrix_, col); it; ++it) {
        out[idx++] = it.row();
    }
}

void BoundaryMatrix::getRowIndices(int row, std::vector<int>& out) const {
    out.clear();
    out.reserve(matrix_.cols());
    if (row < 0 || row >= matrix_.rows()) {
        return;
    }
    for (int j = 0; j < matrix_.cols(); ++j) {
        if (matrix_.coeff(row, j) != 0) {
            out.push_back(j);
        }
    }
    out.shrink_to_fit();
}

const std::vector<int>& BoundaryMatrix::getDimensions() const {
    return dimensions_;
} 