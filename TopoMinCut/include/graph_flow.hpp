#pragma once
#include <vector>

// Computes the min-cut for the given birth/death nodes, generating sets, conflicts, and alphas.
// Returns a pair of vectors: (cut_birth_nodes, cut_death_nodes)
std::pair<std::vector<int>, std::vector<int>>
compute_min_cut(
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
); 

std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_with_dimensions(
    const std::vector<int>& dimensions,
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
); 

std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_simple_cost_function(
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
); 

std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_simple_cost_function_with_dimensions(
    const std::vector<int>& dimensions,
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
);