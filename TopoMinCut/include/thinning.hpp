#pragma once

#include "boundary_matrix.hpp"
#include <queue>
#include <vector>
#include <unordered_map>

class Thinning {
public:
    Thinning(const BoundaryMatrix& boundary_matrix);
    
    // Perform foreground thinning with protected indices (foreground generators)
    void performForegroundThinning(const std::vector<int>& protected_indices);
    
    // Perform background thinning
    void performBackgroundThinning(const std::vector<int>& protected_indices);
    
    // Get the resulting thinned matrix
    const Eigen::SparseMatrix<int>& getThinnedMatrix() const;
    
    // Get the list of removed pairs {child, parent}
    const std::vector<std::pair<int, int>>& getRemovedPairs() const;
    
    // Get remaining negative indices after thinning
    std::vector<int> getRemainingNegativeIndices() const;
    
    // Get remaining positive indices after thinning
    std::vector<int> getRemainingPositiveIndices() const;

private:
    const BoundaryMatrix& boundary_matrix_;
    Eigen::SparseMatrix<int> thinned_matrix_;
    std::vector<std::pair<int, int>> removed_pairs_;
    
    // Helper methods
    std::pair<Eigen::SparseMatrix<int>, std::vector<int>> getParentsMatrix(int max_index) const;
    std::pair<Eigen::SparseMatrix<int>, std::vector<int>> getParentsMatrixBackground(int max_index) const;
    bool isCollapsible(int cell_index) const;
    void collapseCell(int cell_index);
    
    // Priority queue element type
    struct QueueElement {
        double value;
        int child;
        int parent;
        
        // For priority queue ordering (smaller value = higher priority)
        bool operator<(const QueueElement& other) const {
            return value > other.value;  // For min heap (pop lowest first)
        }
    };
    
    // Data structures for tracking relationships
    std::vector<std::vector<int>> row_lists_;  // For each column, list of rows with non-zero entries
    std::vector<std::vector<int>> col_lists_;  // For each row, list of columns with non-zero entries
    std::vector<bool> active_;                 // Track active/inactive status of cells
}; 