#include "iset.hpp"
#include <unordered_set>
#include <algorithm>

std::vector<int> compute_iset_indices(const std::vector<int>& cut_indices, int total_size) {
    std::unordered_set<int> cut_set(cut_indices.begin(), cut_indices.end());
    std::vector<int> iset;
    for (int i = 0; i < total_size; ++i) {
        if (cut_set.find(i) == cut_set.end()) {
            iset.push_back(i);
        }
    }
    return iset;
}

std::vector<int> union_generating_sets(const std::vector<std::vector<int>>& generating_sets, const std::vector<int>& iset_indices) {
    std::unordered_set<int> union_set;
    for (int idx : iset_indices) {
        for (int cell : generating_sets[idx]) {
            union_set.insert(cell);
        }
    }
    return std::vector<int>(union_set.begin(), union_set.end());
} 