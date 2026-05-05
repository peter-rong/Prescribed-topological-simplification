#pragma once

#include "boundary_matrix.hpp"
#include <Eigen/Sparse>
#include <limits>
#include <string>
#include <vector>

namespace topomincut {

struct RunParams {
    int topK = 0;
    double core = -std::numeric_limits<double>::infinity();
    double neighborhood = std::numeric_limits<double>::infinity();
    int cavitySkip = 0;
    int handleSkip = 0;
    int componentSkip = 0;
    /// Safety cap on persistence–min-cut rounds (each round may reorder columns).
    int maxTopologyIterations = 256;
};

/// Populated when non-null pointer is passed to runFromBoundaryMatrix / runFromEigenSparse.
struct RunOutputs {
    /// Alphas in **input** column index order (same indexing as the boundary matrix passed into runFromBoundaryMatrix).
    std::vector<double> alphas_updated;
    std::vector<int> remaining_negative;
    std::vector<int> protected_indices;
    /// Last iteration’s column-permuted boundary matrix (sparse square).
    Eigen::SparseMatrix<int> permuted_boundary_matrix;
    /// True if more than one persistence round ran (legacy “second iteration” path).
    bool completed_second_iteration = false;
    /// Number of persistence / min-cut rounds executed (≥ 1).
    int iterations_completed = 0;
};

/// Full TopoMinCut pipeline in memory only (no disk writes in this translation unit).
int runFromBoundaryMatrix(BoundaryMatrix& boundary_matrix, const RunParams& params, RunOutputs* outs = nullptr);

/// Convenience: packages Eigen + dimensions + alphas then runs the pipeline.
int runFromEigenSparse(const Eigen::SparseMatrix<int>& matrix,
                       const std::vector<int>& dimensions,
                       const std::vector<double>& alphas,
                       const RunParams& params,
                       RunOutputs* outs = nullptr);

} // namespace topomincut
