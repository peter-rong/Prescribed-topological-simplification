#pragma once

#include "boundary_matrix.hpp"
#include <phat/representations/bit_tree_pivot_column.h>
#include <phat/algorithms/exhaustive_compress_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/compute_persistence_pairs.h>
#include <phat/helpers/dualize.h>
#include <array>
#include <vector>

// Structure to hold persistence pair with alpha values and its generator
struct PersistencePair {
    int birth_cell;
    int death_cell;
    double birth_alpha;
    double death_alpha;
    std::vector<int> generator;  // Non-zero entries in death cell's column of reduced matrix
};

class Persistence {
public:
    Persistence(const BoundaryMatrix& boundary_matrix, const std::string& matrix_file);
    ~Persistence();
    
    // Compute persistent homology
    void computePersistence(double core, double neighborhood);
    
    // Get the filtered persistence pairs (with one positive and one negative alpha)
    const std::vector<PersistencePair>& getFilteredPairs() const;
    
    // Get the generator for a given dimension
    const std::array<int, 1000>& getGenerator(int dimension) const;  // Fixed size array for generators
    size_t getGeneratorSize() const;  // Get actual number of generators

    // Get birth and death cell indices
    const std::vector<int>& getBirthCells() const;
    const std::vector<int>& getDeathCells() const;

    // Get the reduced matrices
    const Eigen::SparseMatrix<int>& getReducedMatrix() const;
    const Eigen::SparseMatrix<int>& getDualReducedMatrix() const;
    
    // Get the PHAT matrices
    const phat::boundary_matrix<phat::bit_tree_pivot_column>& getDualPHATMatrix() const;
    
    std::vector<double> outputPersistenceDiagram(const std::string& filename, double minX, double maxX, double minY, double maxY) const;

private:
    const BoundaryMatrix& boundary_matrix_;
    std::string matrix_file_;  // Store the matrix file path
    std::vector<PersistencePair> filtered_pairs_;  // Store filtered persistence pairs
    std::vector<int> birth_cells_;    // Store birth cell indices
    std::vector<int> death_cells_;    // Store death cell indices
    std::array<int, 1000> generators_;  // Fixed size array for generators
    size_t generator_size_;  // Actual number of generators
    Eigen::SparseMatrix<int> reduced_matrix_;
    Eigen::SparseMatrix<int> dual_reduced_matrix_;
    
    // Helper methods
    void loadFromFile(const std::string& filename);
    /// Populate PHAT matrices from BoundaryMatrix::getMatrix() (no disk IO).
    void loadFromBoundaryMatrix(bool verbose_compare_with_eigen);
    void performExhaustiveReduction();
    void performDualReduction();  // Only for getting dual reduced matrix
    void convertFromPHAT();
    void filterPairs(double core, double neighborhood);  // Filter pairs based on alpha values
    
    // PHAT data structures
    phat::boundary_matrix<phat::bit_tree_pivot_column> phat_matrix_;
    phat::boundary_matrix<phat::bit_tree_pivot_column> dual_phat_matrix_;  // For dual reduction
}; 