#include "tags.hpp"
#include <algorithm>
#include <queue>
#include <iostream>
#include <unordered_set>

Tags::Tags(const BoundaryMatrix& boundary_matrix)
    : boundary_matrix_(boundary_matrix) {
    
    // Initialize tags array to 0
    tags_.resize(boundary_matrix_.getAlphaSize(), -1);
    
    // Initialize matrices
    initializeMatrices();
}

void Tags::initializeTags(const std::vector<std::pair<int, int>>& foreground_pairs,
                         const std::vector<std::pair<int, int>>& background_pairs) {
    
    // Fill tags from higher to lower for foreground pairs
    for (const auto& pair : foreground_pairs) {
        int parent = pair.first;
        int child = pair.second;
        tags_[child] = parent;
    }
    
    // Fill tags from higher to lower for background pairs
    for (const auto& pair : background_pairs) {
        int parent = pair.first;
        int child = pair.second;
        tags_[child] = parent;
    }
}

void Tags::initializeMatrices() {
    
    int alpha_size = boundary_matrix_.getAlphaSize();
    const auto& matrix = boundary_matrix_.getMatrix();
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    
    // Find the first non-negative alpha index
    int first_non_negative = 0;
    while (first_non_negative < alpha_size && alphas[first_non_negative] <= 0) {
        first_non_negative++;
    }
    std::cout << "first_non_negative: " << first_non_negative << std::endl;
    
    // Initialize matrices
    fore_matrix_ = matrix;
    back_matrix_ = matrix;
    
    // For foreground matrix: clear rows and columns >= first_non_negative
    for (int k = 0; k < fore_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(fore_matrix_, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (row >= first_non_negative || col >= first_non_negative) {
                it.valueRef() = 0;
            }
        }
    }
    fore_matrix_.prune([](int, int, int value) { return value != 0; });  // Remove all zero elements
    
    // For background matrix: clear rows and columns < first_non_negative
    // We want to keep elements where BOTH row and col are >= first_non_negative
    for (int k = 0; k < back_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(back_matrix_, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            if (row < first_non_negative || col < first_non_negative) {
                it.valueRef() = 0;
            }
        }
    }
    back_matrix_.prune([](int, int, int value) { return value != 0; });  // Remove all zero elements
    
    // Make matrices compressed
    fore_matrix_.makeCompressed();
    back_matrix_.makeCompressed();
    
    // Print matrix statistics
    std::cout << "\nForeground matrix statistics:" << std::endl;
    std::cout << "  Size: " << fore_matrix_.rows() << "x" << fore_matrix_.cols() << std::endl;
    std::cout << "  Non-zeros: " << fore_matrix_.nonZeros() << std::endl;
    std::cout << "  Density: " << (double)fore_matrix_.nonZeros() / (fore_matrix_.rows() * fore_matrix_.cols()) * 100 << "%" << std::endl;
    
    std::cout << "\nBackground matrix statistics:" << std::endl;
    std::cout << "  Size: " << back_matrix_.rows() << "x" << back_matrix_.cols() << std::endl;
    std::cout << "  Non-zeros: " << back_matrix_.nonZeros() << std::endl;
    std::cout << "  Density: " << (double)back_matrix_.nonZeros() / (back_matrix_.rows() * back_matrix_.cols()) * 100 << "%" << std::endl;
}

void Tags::initializeBirthDeathInfo(const std::vector<std::pair<int, int>>& filtered_pairs) {
    
    // Get alpha values
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    
    // Reserve space for all pairs
    birth_cells_.reserve(filtered_pairs.size());
    death_cells_.reserve(filtered_pairs.size());
    birth_times_.reserve(filtered_pairs.size());
    death_times_.reserve(filtered_pairs.size());
    
    // Process each filtered pair
    for (const auto& pair : filtered_pairs) {
        int birth = pair.first;
        int death = pair.second;
        double birth_time = std::abs(alphas[birth]);
        double death_time = std::abs(alphas[death]);
        
        birth_cells_.push_back(birth);
        death_cells_.push_back(death);
        birth_times_.push_back(birth_time);
        death_times_.push_back(death_time);
    }
}

std::vector<int> Tags::getGeneratingSetForeFast(int starting_index) const {
    
    // Pre-allocate vectors with reasonable sizes
    std::vector<int> set;
    set.reserve(fore_matrix_.rows());  // Reserve max possible size
    
    std::vector<bool> visited(fore_matrix_.rows(), false);
    std::queue<int> children_queue;
    
    // Initialize with starting index
    children_queue.push(starting_index);
    visited[starting_index] = true;
    
    // Pre-build row and column lists for faster access
    std::vector<std::vector<int>> row_lists(fore_matrix_.rows());
    std::vector<std::vector<int>> col_lists(fore_matrix_.cols());
    
    // Pre-allocate space in lists
    for (auto& list : row_lists) list.reserve(100);  // Reasonable initial capacity
    for (auto& list : col_lists) list.reserve(100);
    
    // Build lists in a single pass
    for (int k = 0; k < fore_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(fore_matrix_, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            row_lists[col].push_back(row);
            col_lists[row].push_back(col);
        }
    }
    
    // Pre-allocate temporary vectors for parents and children
    std::vector<int> curr_parents;
    curr_parents.reserve(100);
    std::vector<int> next_children;
    next_children.reserve(100);
    
    while (!children_queue.empty()) {
        int curr_child = children_queue.front();
        children_queue.pop();
        
        set.push_back(curr_child);
        
        // Get unvisited parents
        curr_parents.clear();
        for (int parent : col_lists[curr_child]) {
            if (!visited[parent]) {
                curr_parents.push_back(parent);
                visited[parent] = true;
            }
        }
        
        // Add unvisited parents to queue
        for (int parent : curr_parents) {
            children_queue.push(parent);
        }
        
        // Get next children from tags
        next_children.clear();
        for (int parent : curr_parents) {
            int child = tags_[parent];
            if (child != -1 && !visited[child]) {  // Check if tag exists and not visited
                next_children.push_back(child);
                visited[child] = true;
            }
        }
        
        // Add unvisited next children to queue
        for (int child : next_children) {
            children_queue.push(child);
        }
    }
    
    return set;
}

std::vector<int> Tags::getGeneratingSetBackFast(int starting_index) const {
    
    std::vector<std::vector<int>> row_lists(back_matrix_.rows());
    std::vector<std::vector<int>> col_lists(back_matrix_.cols());
    
    // Pre-allocate space in lists
    for (auto& list : row_lists) list.reserve(100);  // Reasonable initial capacity
    for (auto& list : col_lists) list.reserve(100);
    
    // Build lists in a single pass
    for (int k = 0; k < back_matrix_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(back_matrix_, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            row_lists[col].push_back(row);
            col_lists[row].push_back(col);
        }
    }
    
    // Pre-allocate result vector with maximum possible size
    std::vector<int> result;
    result.reserve(back_matrix_.rows());
    
    // Pre-allocate temporary vectors for parents and children
    std::vector<int> curr_parents;
    curr_parents.reserve(100);
    std::vector<int> next_children;
    next_children.reserve(100);
    
    std::vector<bool> visited(back_matrix_.rows(), false);
    std::queue<int> children_queue;
    
    // Initialize with starting index
    children_queue.push(starting_index);
    visited[starting_index] = true;
    
    while (!children_queue.empty()) {
        int curr_child = children_queue.front();
        children_queue.pop();
        
        result.push_back(curr_child);
        
        // Get unvisited parents
        curr_parents.clear();
        for (int parent : row_lists[curr_child]) {
            if (!visited[parent]) {
                curr_parents.push_back(parent);
                visited[parent] = true;
            }
        }
        
        // Add unvisited parents to queue
        for (int parent : curr_parents) {
            children_queue.push(parent);
        }
        
        // Get next children from tags
        next_children.clear();
        for (int parent : curr_parents) {
            int child = tags_[parent];
            if (child != -1 && !visited[child]) {  // Check if tag exists and not visited
                next_children.push_back(child);
                visited[child] = true;
            }
        }
        
        // Add unvisited next children to queue
        for (int child : next_children) {
            children_queue.push(child);
        }
    }

    return result;
}

std::vector<int> Tags::getGeneratingSetForeFast2(const std::vector<std::vector<int>>& col_lists,
                                               int starting_index) const {
    
    // Pre-allocate vectors with reasonable sizes
    std::vector<int> set;
    set.reserve(fore_matrix_.rows());  // Reserve max possible size
    
    std::vector<bool> visited(fore_matrix_.rows(), false);
    std::queue<int> children_queue;
    
    // Initialize with starting index
    children_queue.push(starting_index);
    visited[starting_index] = true;
    
    while (!children_queue.empty()) {
        int curr_child = children_queue.front();
        children_queue.pop();
        
        set.push_back(curr_child);
        
        // Process parents and their tagged children in a single pass
        for (int parent : col_lists[curr_child]) {
            if (!visited[parent]) {
                visited[parent] = true;
                children_queue.push(parent);
                
                // Immediately process tagged child if it exists and is unvisited
                int tagged_child = tags_[parent];
                if (tagged_child != -1 && !visited[tagged_child]) {
                    visited[tagged_child] = true;
                    children_queue.push(tagged_child);
                }
            }
        }
    }
    
    return set;
}

std::vector<int> Tags::getGeneratingSetBackFast2(const std::vector<std::vector<int>>& row_lists,
                                               int starting_index) const {
    
    // Pre-allocate result vector with maximum possible size
    std::vector<int> result;
    result.reserve(back_matrix_.rows());
    
    std::vector<bool> visited(back_matrix_.rows(), false);
    std::queue<int> children_queue;
    
    // Initialize with starting index
    children_queue.push(starting_index);
    visited[starting_index] = true;
    
    while (!children_queue.empty()) {
        int curr_child = children_queue.front();
        children_queue.pop();
        
        result.push_back(curr_child);
        
        // Process parents and their tagged children in a single pass
        for (int parent : row_lists[curr_child]) {
            if (!visited[parent]) {
                visited[parent] = true;
                children_queue.push(parent);
                
                // Immediately process tagged child if it exists and is unvisited
                int tagged_child = tags_[parent];
                if (tagged_child != -1 && !visited[tagged_child]) {
                    visited[tagged_child] = true;
                    children_queue.push(tagged_child);
                }
            }
        }
    }

    return result;
}

std::vector<std::vector<int>> Tags::computeConflicts(
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices) const {
    
    // Initialize conflict matrix with zeros
    std::vector<std::vector<int>> conflicts(real_birth_cells.size(), 
                                          std::vector<int>(real_death_cells.size(), 0));
    
    // Create foreground tags array - each cell records which generating sets it's in
    int alpha_size = boundary_matrix_.getAlphaSize();
    
    // Use vector of vectors instead of raw arrays
    std::vector<std::vector<int>> fore_tags(alpha_size);
    
    // Process birth cells
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        const auto& generating_sets_fore_indices = real_generating_sets_fore_indices[i];
        for (int idx : generating_sets_fore_indices) {
            // Record that this cell is in generating set i
            fore_tags[idx].push_back(i);
        }
    }
    
    // Process each death cell
    for (size_t i = 0; i < real_death_cells.size(); ++i) {
        const auto& generating_sets_back_indices = real_generating_sets_back_indices[i];
        std::unordered_set<int> conflicting_fore_sets;
        
        // Process each cell in the generating set
        for (int curr : generating_sets_back_indices) {
            // Get neighbors using Eigen's sparse matrix iterator
            for (Eigen::SparseMatrix<int>::InnerIterator it(boundary_matrix_.getMatrix(), curr); it; ++it) {
                int neighbor = it.row();
                // Add all foreground sets this neighbor is in to conflicting sets
                for (int fore_set : fore_tags[neighbor]) {
                    conflicting_fore_sets.insert(fore_set);
                }
            }
            
            // Also check tagged child if exists
            int child = tags_[curr];
            if (child != -1) {
                // Add all foreground sets this child is in to conflicting sets
                for (int fore_set : fore_tags[child]) {
                    conflicting_fore_sets.insert(fore_set);
                }
            }
        }
        
        // Mark conflicts for this death cell
        for (int fore_set : conflicting_fore_sets) {
            conflicts[fore_set][i] = 1;
        }
    }
    
    return conflicts;
} 
