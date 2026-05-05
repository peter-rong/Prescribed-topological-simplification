#include "thinning.hpp"
#include <algorithm>
#include <queue>
#include <iostream>
#include <set>

Thinning::Thinning(const BoundaryMatrix& boundary_matrix)
    : boundary_matrix_(boundary_matrix) {
    // Initialize thinned matrix as a copy of the boundary matrix
    thinned_matrix_ = boundary_matrix.getMatrix();
}

std::pair<Eigen::SparseMatrix<int>, std::vector<int>> 
Thinning::getParentsMatrix(int max_index) const {
    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Create a new sparse matrix with the same size
    Eigen::SparseMatrix<int> sub_matrix(alpha_size, alpha_size);
    sub_matrix.reserve(thinned_matrix_.nonZeros());
    
    // Copy only the foreground region (before max_index)
    for (int k = 0; k < thinned_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(thinned_matrix_, k); it; ++it) {
            if (it.row() < max_index && it.col() < max_index) {
                sub_matrix.insert(it.row(), it.col()) = it.value();
            }
        }
    }
    sub_matrix.makeCompressed();
    
    // Pre-allocate and count non-zeros in each column
    std::vector<int> counts(alpha_size, 0);
    for (int k = 0; k < sub_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(sub_matrix, k); it; ++it) {
            counts[it.row()]++;
        }
    }

    return {sub_matrix, counts};
}

std::pair<Eigen::SparseMatrix<int>, std::vector<int>> 
Thinning::getParentsMatrixBackground(int max_index) const {

    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Create a new sparse matrix with the same size
    Eigen::SparseMatrix<int> sub_matrix(alpha_size, alpha_size);
    sub_matrix.reserve(thinned_matrix_.nonZeros());
    
    // Copy only the background region (after max_index)
    for (int k = 0; k < thinned_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(thinned_matrix_, k); it; ++it) {
            if (it.row() >= max_index && it.col() >= max_index) {
                sub_matrix.insert(it.row(), it.col()) = it.value();
            }
        }
    }
    sub_matrix.makeCompressed();
    
    // Pre-allocate and count non-zeros in each column
    std::vector<int> counts(alpha_size, 0);
    for (int k = 0; k < sub_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(sub_matrix, k); it; ++it) {
            counts[it.col()]++;
        }
    }

    return {sub_matrix, counts};
}

void Thinning::performForegroundThinning(const std::vector<int>& protected_indices) {
    // Store protected indices
    // Get alpha values
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Find the maximum index where alpha value is non-positive
    int max_index = 0;
    for (int i = 0; i < alpha_size; ++i) {
        if (alphas[i] <= 0) {
            max_index = i + 1;
        }
    }
    
    std::cout << "Max negative index: " << max_index << std::endl;
    
    // Get parent matrix and counts
    auto result = getParentsMatrix(max_index);
    const auto& temp_matrix = result.first;
    std::vector<int> counts = result.second;
    
    // Pre-allocate all data structures with exact sizes
    active_.resize(max_index, true);
    row_lists_.resize(max_index);
    col_lists_.resize(max_index);
    removed_pairs_.clear();
    removed_pairs_.reserve(max_index); // Reserve space for worst case
    
    // Set protected indices as inactive
    for (int idx : protected_indices) {
        active_[idx] = false;
    }

    for (size_t i = max_index; i < active_.size(); ++i) {
        active_[i] = false;
    }
    // Count non-zeros for each row and column
    std::vector<int> row_counts(max_index, 0);
    std::vector<int> col_counts(max_index, 0);
    
    // First pass: count non-zeros
    for (int k = 0; k < temp_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(temp_matrix, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (row < max_index && col < max_index) {
                row_counts[col]++;
                col_counts[row]++;
            }
        }
    }
    
    // Pre-allocate row and column lists with exact sizes
    for (int i = 0; i < max_index; ++i) {
        row_lists_[i].resize(row_counts[i]);
        col_lists_[i].resize(col_counts[i]);
        row_counts[i] = 0;  // Reset for second pass
        col_counts[i] = 0;  // Reset for second pass
    }
    
    // Second pass: populate lists
    for (int k = 0; k < temp_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(temp_matrix, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (row < max_index && col < max_index) {
                row_lists_[col][row_counts[col]++] = row;
                col_lists_[row][col_counts[row]++] = col;
            }
        }
    }
    
    // Create priority queue
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::less<QueueElement>> heap;
    
    // Pre-allocate vectors for children and parents
    std::vector<int> children(max_index);
    std::vector<int> active_children(max_index);
    size_t num_children = 0;
    size_t num_active_children = 0;
    
    // Populate heap with initial simple pairs
    for (int i = 0; i < max_index; ++i) {
        if (active_[i] && counts[i] == 0) {
            // Find children with count 1
            num_children = 0;
            for (int row : row_lists_[i]) {
                if (counts[row] == 1) {
                    children[num_children++] = row;
                }
            }
            
            // Filter active children and check if they are isolated
            num_active_children = 0;
            for (size_t j = 0; j < num_children; ++j) {
                int child = children[j];
                if (active_[child]) {
                    // Check if child is isolated (has no other active parents)
                    bool is_isolated = true;
                    for (int p : col_lists_[child]) {
                        if (p != i && active_[p]) {
                            is_isolated = false;
                            break;
                        }
                    }
                    if (is_isolated) {
                        active_children[num_active_children++] = child;
                    }
                }
            }
            
            if (num_active_children > 0) {

                int child = active_children[num_active_children-1];
                
                // Don't mark as inactive yet - will be marked when actually removed
                heap.push({static_cast<double>(-1*alphas[i]), child, i});
                //heap.push({std::abs(alphas[i]), child, i});
            }
        }
    }
    
    // Pre-allocate vector for parents
    std::vector<int> parents(max_index);
    size_t num_parents = 0;
    
    // Main thinning loop: process one pair at a time
    while (!heap.empty()) {

        double current_value = heap.top().value;

        // Collect all pairs with the same minimal priority
        std::vector<QueueElement> batch;
        while (!heap.empty() && heap.top().value == current_value) {
            batch.push_back(heap.top());
            heap.pop();
        }
        for (const auto& selected : batch) {

            int child = selected.child;
            int parent = selected.parent;
            
            // Verify this pair is still a valid simple pair before removing
            // Check: parent is active, has count 0, child is active, has count 1, and is isolated
            if (!active_[parent] || !active_[child] || counts[parent] != 0 || counts[child] != 1) {
                continue; // Skip invalid pairs
            }
            
            // Check if child is still isolated (has no other active parents)
            bool is_isolated = true;
            for (int p : col_lists_[child]) {
                if (p != parent && active_[p]) {
                    is_isolated = false;
                    break;
                }
            }
            if (!is_isolated) {
                continue; // Skip if child is no longer isolated
            }
            
            // Now it's safe to remove this pair
            active_[child] = false;
            active_[parent] = false;
            removed_pairs_.emplace_back(child, parent);
            
            // Update children of removed cells
            const auto& children1 = row_lists_[child];
            const auto& children2 = row_lists_[parent];
            
            // Update counts in a single pass
            for (int c : children1) counts[c]--;
            for (int c : children2) counts[c]--;
            
            // Process children function with pre-allocated vectors
            auto processChildren = [&](const auto& children) {
                for (int c : children) {
                    if (active_[c] && counts[c] == 1) {
                        // Find active parents
                        num_parents = 0;
                        for (int p : col_lists_[c]) {
                            if (active_[p]) {
                                parents[num_parents++] = p;
                            }
                        }
                        
                        if (num_parents == 1) {
                            int parent_local = parents[0];
                            if (active_[parent_local] && counts[parent_local] == 0) {
                                // Don't mark as inactive yet - will be marked when actually removed
                                heap.push({static_cast<double>(-1 * alphas[parent_local]), c, parent_local});
                            }
                        }
                    }
                }
            };
            
            processChildren(children1);
            processChildren(children2);
        }
    }
    
}

void Thinning::performBackgroundThinning(const std::vector<int>& protected_indices) {
    // Get alpha values
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Find the minimum index where alpha value is positive
    int min_index = alpha_size;
    for (int i = 0; i < alpha_size; ++i) {
        if (alphas[i] > 0) {
            min_index = i;
            break;
        }
    }
    
    // Get parent matrix and counts for positive indices
    auto result = getParentsMatrixBackground(min_index);
    const auto& temp_matrix = result.first;
    std::vector<int> counts = result.second;
    
    // Pre-allocate all data structures with exact sizes
    active_.resize(alpha_size, true);
    row_lists_.resize(alpha_size);
    col_lists_.resize(alpha_size);
    removed_pairs_.clear();
    removed_pairs_.reserve(alpha_size);
    
    // Set protected indices as inactive
    for (int idx : protected_indices) {
        active_[idx] = false;
    }

    for (size_t i = 0; i < min_index; ++i) {
        active_[i] = false;
    }
    
    // Count non-zeros for each row and column
    std::vector<int> row_counts(alpha_size, 0);
    std::vector<int> col_counts(alpha_size, 0);
    
    // First pass: count non-zeros
    for (int k = 0; k < temp_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(temp_matrix, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            row_counts[col]++;
            col_counts[row]++;
        }
    }
    
    // Pre-allocate row and column lists with exact sizes
    for (int i = min_index; i < alpha_size; ++i) {
        row_lists_[i].resize(row_counts[i]);
        col_lists_[i].resize(col_counts[i]);
        row_counts[i] = 0;  // Reset for second pass
        col_counts[i] = 0;  // Reset for second pass
    }
    
    // Second pass: populate lists
    for (int k = 0; k < temp_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(temp_matrix, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            row_lists_[col][row_counts[col]++] = row;
            col_lists_[row][col_counts[row]++] = col;
        }
    }
    
    // Create priority queue
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::less<QueueElement>> heap;
    
    // Pre-allocate vectors for children and parents
    std::vector<int> children(alpha_size);
    std::vector<int> active_children(alpha_size);
    size_t num_children = 0;
    size_t num_active_children = 0;
    
    // Populate heap with initial simple pairs
    for (int i = 0; i < alpha_size; ++i) {
        if (active_[i] && counts[i] == 0) {
            // Find children with count 1
            num_children = 0;
            for (int col : col_lists_[i]) {
                if (counts[col] == 1) {
                    children[num_children++] = col;
                }
            }
            
            // Filter active children
            num_active_children = 0;
            for (size_t j = 0; j < num_children; ++j) {
                int child = children[j];
                if (active_[child]) {
                    active_children[num_active_children++] = child;
                }
            }
            
            if (num_active_children > 0) {
                // Find child with smallest alpha value
                int child = active_children[0];
                
                // Don't mark as inactive yet - will be marked when actually removed
                heap.push({alphas[i], child, i});
            }
        }
    }
    
    // Pre-allocate vector for parents
    std::vector<int> parents(alpha_size);
    size_t num_parents = 0;
    
    // Main thinning loop: process one pair at a time
    while (!heap.empty()) {

        double current_value = heap.top().value;

        // Collect all pairs with the same minimal priority
        std::vector<QueueElement> batch;
        while (!heap.empty() && heap.top().value == current_value) {
            batch.push_back(heap.top());
            heap.pop();
        }

        for (const auto& selected : batch) {
            
            int child = selected.child;
            int parent = selected.parent;
            
            // Verify this pair is still a valid simple pair before removing
            // Check: parent is active, has count 0, child is active, has count 1
            if (!active_[parent] || !active_[child] || counts[parent] != 0 || counts[child] != 1) {
                continue; // Skip invalid pairs
            }
            
            // Now it's safe to remove this pair
            active_[child] = false;
            active_[parent] = false;
            removed_pairs_.emplace_back(child, parent);
            
            // Update children of removed cell
            const auto& children1 = col_lists_[child];
            const auto& children2 = col_lists_[parent];
            
            // Update counts in a single pass
            for (int c : children1) counts[c]--;
            for (int c : children2) counts[c]--;
            
            // Process children function with pre-allocated vectors
            auto processChildren = [&](const auto& children) {
                for (int c : children) {
                    if (active_[c] && counts[c] == 1) {
                        // Find active parents
                        num_parents = 0;
                        for (int p : row_lists_[c]) {
                            if (active_[p]) {
                                parents[num_parents++] = p;
                            }
                        }
                        
                        if (num_parents == 1) {
                            int parent_local = parents[0];
                            if (active_[parent_local] && counts[parent_local] == 0) {
                                // Don't mark as inactive yet - will be marked when actually removed
                                heap.push({alphas[parent_local], c, parent_local});
                            }
                        }
                    }
                }
            };
            
            processChildren(children1);
            processChildren(children2);
        }
    }
    
    // Count remaining positive indices
    int remaining_positive = 0;
    for (int i = min_index; i < alpha_size; ++i) {
        if (active_[i]) {
            remaining_positive++;
        }
    }
}

const Eigen::SparseMatrix<int>& Thinning::getThinnedMatrix() const {
    return thinned_matrix_;
}

const std::vector<std::pair<int, int>>& Thinning::getRemovedPairs() const {
    return removed_pairs_;
}

bool Thinning::isCollapsible(int cell_index) const {
    // TODO: Implement collapsibility check
    return false;
}

void Thinning::collapseCell(int cell_index) {
    // TODO: Implement cell collapse
}
 
std::vector<int> Thinning::getRemainingNegativeIndices()const {
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Find max_index (largest index where alpha ≤ 0)
    int max_index = 0;
    for (int i = 0; i < alpha_size; ++i) {
        if (alphas[i] <= 0) {
            max_index = i + 1;
        }
    }
    
    // Pre-allocate result vector with worst-case size
    std::vector<int> remaining;
    remaining.reserve(max_index);
    
    // Create a set of removed indices for faster lookup
    std::set<int> removed_indices;
    for (const auto& pair : removed_pairs_) {
        removed_indices.insert(pair.first);  // child
        removed_indices.insert(pair.second); // parent
    }
    
    // Get all indices up to max_index that are not in removed_pairs
    for (int i = 0; i < max_index; ++i) {
        if (removed_indices.find(i) == removed_indices.end()) {
            remaining.push_back(i);
        }
    }
    
    return remaining;
}

std::vector<int> Thinning::getRemainingPositiveIndices() const {
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    int alpha_size = boundary_matrix_.getAlphaSize();
    // Find min_index (smallest index where alpha > 0)
    int min_index = alpha_size;
    for (int i = 0; i < alpha_size; ++i) {
        if (alphas[i] > 0) {
            min_index = i;
            break;
        }
    }
    // Pre-allocate result vector with worst-case size
    std::vector<int> remaining;
    remaining.reserve(alpha_size - min_index);
    // Create a set of removed indices for faster lookup
    std::set<int> removed_indices;
    for (const auto& pair : removed_pairs_) {
        removed_indices.insert(pair.first);  // child
        removed_indices.insert(pair.second); // parent
    }
    // Get all indices from min_index to alpha_size that are not in removed_pairs
    for (int i = min_index; i < alpha_size; ++i) {
        if (removed_indices.find(i) == removed_indices.end()) {
            remaining.push_back(i);
        }
    }
    return remaining;
}