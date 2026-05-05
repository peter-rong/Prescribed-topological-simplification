#pragma once
#include <vector>

// Computes the complement indices (ISET) given the cut indices and total size
std::vector<int> compute_iset_indices(const std::vector<int>& cut_indices, int total_size);

// Unions the generating sets for the given ISET indices
std::vector<int> union_generating_sets(const std::vector<std::vector<int>>& generating_sets, const std::vector<int>& iset_indices); 