#pragma once

#include "boundary_matrix.hpp"
#include <vector>
#include <string>
#include <Eigen/Sparse>

class Tags {
public:
    Tags(const BoundaryMatrix& boundary_matrix);
    
    // Initialize tags based on removed pairs from foreground and background thinning
    void initializeTags(const std::vector<std::pair<int, int>>& foreground_pairs,
                       const std::vector<std::pair<int, int>>& background_pairs);
    
    // Initialize birth/death info from filtered pairs
    void initializeBirthDeathInfo(const std::vector<std::pair<int, int>>& filtered_pairs);
    
    // Get generating set for a birth cell in foreground
    std::vector<int> getGeneratingSetForeFast(int starting_index) const;
    
    // Get generating set for a birth cell in background
    std::vector<int> getGeneratingSetBackFast(int starting_index) const;
    
    // Get generating set for a birth cell in foreground (with pre-computed lists)
    std::vector<int> getGeneratingSetForeFast2(const std::vector<std::vector<int>>& row_lists, 
                                             int starting_index) const;
    
    // Get generating set for a birth cell in background (with pre-computed lists)
    std::vector<int> getGeneratingSetBackFast2(const std::vector<std::vector<int>>& col_lists,
                                             int starting_index) const;
    
    // Get the tags array
    const std::vector<int>& getTags() const { return tags_; }
    
    // Get the foreground matrix (with non-negative indices zeroed out)
    const Eigen::SparseMatrix<int>& getForegroundMatrix() const { return fore_matrix_; }
    
    // Get the background matrix (with negative indices zeroed out)
    const Eigen::SparseMatrix<int>& getBackgroundMatrix() const { return back_matrix_; }
    
    // Get birth cells, death cells, and their times
    const std::vector<int>& getBirthCells() const { return birth_cells_; }
    const std::vector<int>& getDeathCells() const { return death_cells_; }
    const std::vector<double>& getBirthTimes() const { return birth_times_; }
    const std::vector<double>& getDeathTimes() const { return death_times_; }
    
    // New function to compute conflicts
    
    
    std::vector<std::vector<int>> computeConflicts(
        const std::vector<int>& real_birth_cells,
        const std::vector<int>& real_death_cells,
        const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
        const std::vector<std::vector<int>>& real_generating_sets_back_indices) const;
    
private:
    BoundaryMatrix boundary_matrix_;
    std::vector<int> tags_;  // Tags array initialized to 0
    
    // Matrices for isolation checks
    Eigen::SparseMatrix<int> fore_matrix_;  // Foreground matrix with non-negative indices zeroed
    Eigen::SparseMatrix<int> back_matrix_;  // Background matrix with negative indices zeroed
    
    // Birth and death information
    std::vector<int> birth_cells_;
    std::vector<int> death_cells_;
    std::vector<double> birth_times_;
    std::vector<double> death_times_;
    
    // Helper methods
    void initializeMatrices();
};