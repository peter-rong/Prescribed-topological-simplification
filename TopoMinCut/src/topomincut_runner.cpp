#include "topomincut_runner.hpp"
#include "boundary_matrix.hpp"
#include "persistence.hpp"
#include "thinning.hpp"
#include "tags.hpp"
#include "iset.hpp"
#include "graph_flow.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <memory>
#include <limits>
#include <queue>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>
#include <unordered_set>

namespace topomincut {

int runFromBoundaryMatrix(BoundaryMatrix& boundary_matrix, const RunParams& params, RunOutputs* outs) {
    std::cout << "Starting TopoMinCut..." << std::endl;

    using namespace std;

    const int topK = params.topK;
    double core = params.core;
    double neighborhood = params.neighborhood;
    const int cavitySkip = params.cavitySkip;
    const int handleSkip = params.handleSkip;
    const int componentSkip = params.componentSkip;

    std::cout << "Matrix dimensions: " << boundary_matrix.getMatrix().rows()
              << " x " << boundary_matrix.getMatrix().cols() << std::endl;
    std::cout << "Number of alpha values: " << boundary_matrix.getAlphaSize() << std::endl;

    const int n0 = boundary_matrix.getAlphaSize();
    std::vector<int> col_original_id(static_cast<size_t>(n0));
    std::iota(col_original_id.begin(), col_original_id.end(), 0);
    const int max_topo_iters = std::max(1, params.maxTopologyIterations);

    for (int iteration = 0;; ++iteration) {
        std::cout << "\n=== TopoMinCut iteration " << (iteration + 1) << " ===" << std::endl;

    Persistence persistence(boundary_matrix, "");
    persistence.computePersistence(core, neighborhood);

    // Get filtered pairs (with one positive and one negative alpha)
    const auto& filtered_pairs = persistence.getFilteredPairs();
    std::cout << "Found " << filtered_pairs.size() << " filtered persistence pairs:\n\n";

    // Sort pairs by survival time (highest to lowest)
    std::vector<PersistencePair> sorted_pairs = filtered_pairs;
    std::sort(sorted_pairs.begin(), sorted_pairs.end(),
        [](const PersistencePair& a, const PersistencePair& b) {
            double survival_a = std::abs(a.death_alpha - a.birth_alpha);
            double survival_b = std::abs(b.death_alpha - b.birth_alpha);
            return survival_a > survival_b;  // Sort in descending order
        });

    // Use a set to store unique standard generators
    std::set<int> protected_indices;
    std::set<int> protected_indices_background;

    // Display all pairs but only use selected ones for generating sets
    for (const auto& pair : sorted_pairs) {
        
        // Print standard generator and add to protected indices
        //std::cout << "Standard generator: [";
        for (size_t i = 0; i < pair.generator.size(); ++i) {
            //std::cout << pair.generator[i];
            //if (i < pair.generator.size() - 1) std::cout << ", ";
            // Add to protected indices
            protected_indices.insert(pair.generator[i]);
        }
        //std::cout << "]\n";
        
        // Get and print dual generator from PHAT matrix using birth cell's row
        std::vector<int> dual_generator;
        const auto& dual_matrix = persistence.getDualPHATMatrix();
        int transformed_birth_cell = boundary_matrix.getAlphaSize() - 1 - pair.birth_cell;
        
        // Get the row data for the transformed birth cell
        phat::column row_data;
        dual_matrix.get_col(transformed_birth_cell, row_data);
        dual_generator = std::vector<int>(row_data.begin(), row_data.end());
        
       //std::cout << "Dual generator: [";
        //for (size_t i = 0; i < dual_generator.size(); ++i) {
            //std::cout << boundary_matrix.getAlphaSize() - 1 - dual_generator[i];
            //if (i < dual_generator.size() - 1) std::cout << ", ";
        //}
        //std::cout << "]\n\n";
        // Add dual generator indices to protected_indices_background
        for (int idx : dual_generator) {
            protected_indices_background.insert(boundary_matrix.getAlphaSize() - 1 - idx);
        }
    }
    
    // Find the last cell with infinity alpha value
    int last_infinity_idx = -1;

    auto alphas = boundary_matrix.getAlphas();

    for (int i = alphas.size() - 1; i >= 0; --i) {
        if (alphas[i] == std::numeric_limits<double>::infinity()) {
            last_infinity_idx = i;
            break;
        }
    }
    if (last_infinity_idx != -1) {
        protected_indices_background.insert(last_infinity_idx);
    }


    //also for every index that has alpha value < core, add to protected_indices
    
    for (int i = 0; i < alphas.size(); ++i) {
        if (alphas[i] < core) {
            protected_indices.insert(i);
        }
    }

    //for every index that has alpha value > neighborhood, add to protected_indices_background
    for (int i = 0; i < alphas.size(); ++i) {
        if (alphas[i] > neighborhood) {
            protected_indices_background.insert(i);
        }
    }
    //now add the last infinity index to the protected indices
    
    std::vector<int> protected_indices_background_vec(protected_indices_background.begin(), protected_indices_background.end());
    
    // Convert set to vector for thinning
    std::vector<int> protected_indices_vec(protected_indices.begin(), protected_indices.end());
    
    std::cout << "Number of protected indices: " << protected_indices_vec.size() << std::endl;

    std::cout << std::endl;

    // Perform foreground thinning
    std::cout << "\nPerforming foreground thinning..." << std::endl;
    Thinning foreground_thinning(boundary_matrix);
    foreground_thinning.performForegroundThinning(protected_indices_vec);
    
    // Print foreground thinning results
    const auto& foreground_removed_pairs = foreground_thinning.getRemovedPairs();

    if (!foreground_removed_pairs.empty()) {

        std::cout << "\nForeground thinning results:" << std::endl;
        std::cout << "Number of removed pairs: " << foreground_removed_pairs.size() << std::endl;
        std::cout << "First thinned pair - Child: " << foreground_removed_pairs.front().first 
                  << " (alpha=" << alphas[foreground_removed_pairs.front().first] 
                  << "), Parent: " << foreground_removed_pairs.front().second 
                  << " (alpha=" << alphas[foreground_removed_pairs.front().second] << ")" << std::endl;
        
        std::cout << "Last thinned pair - Child: " << foreground_removed_pairs.back().first 
                  << " (alpha=" << alphas[foreground_removed_pairs.back().first] 
                  << "), Parent: " << foreground_removed_pairs.back().second 
                  << " (alpha=" << alphas[foreground_removed_pairs.back().second] << ")" << std::endl;
        
    }
    
    // Get remaining negative indices after foreground thinning
    std::vector<int> remaining_negative = foreground_thinning.getRemainingNegativeIndices();
    std::cout << "Number of remaining negative indices after foreground thinning: " << remaining_negative.size() << std::endl;

    
    // Perform background thinning
    std::cout << "\nPerforming background thinning..." << std::endl;
    Thinning background_thinning(boundary_matrix);
    background_thinning.performBackgroundThinning(protected_indices_background_vec);
    
    // Get remaining positive indices after background thinning
    std::vector<int> remaining_positive = background_thinning.getRemainingPositiveIndices();
    std::cout << "Number of remaining positive indices after background thinning: " << remaining_positive.size() << std::endl;
    
    // Print background thinning results
    const auto& background_removed_pairs = background_thinning.getRemovedPairs();
    if (!background_removed_pairs.empty()) {
        
        std::cout << "\nBackground thinning results:" << std::endl;
        std::cout << "Number of removed pairs: " << background_removed_pairs.size() << std::endl;
        std::cout << std::fixed << std::setprecision(7);  // Set fixed precision
        std::cout << "First thinned pair - Child: " << background_removed_pairs.front().first 
                  << " (alpha=" << alphas[background_removed_pairs.front().first] 
                  << "), Parent: " << background_removed_pairs.front().second 
                  << " (alpha=" << alphas[background_removed_pairs.front().second] << ")" << std::endl;
        
        std::cout << "Last thinned pair - Child: " << background_removed_pairs.back().first 
                  << " (alpha=" << alphas[background_removed_pairs.back().first] 
                  << "), Parent: " << background_removed_pairs.back().second 
                  << " (alpha=" << alphas[background_removed_pairs.back().second] << ")" << std::endl;
        
    }
    

    // Create Tags object
    Tags tags(boundary_matrix);
    
    // Initialize tags with removed pairs
    tags.initializeTags(foreground_removed_pairs, background_removed_pairs);
    
    // Convert filtered persistence pairs to the format needed by Tags
    std::vector<std::pair<int, int>> converted_pairs;
    converted_pairs.reserve(sorted_pairs.size());
    for (const auto& pair : sorted_pairs) {
        converted_pairs.emplace_back(pair.birth_cell, pair.death_cell);
    }

    // Initialize birth/death info
    tags.initializeBirthDeathInfo(converted_pairs);
    
    // Get generating sets for all birth cells
    const auto& birth_cells = tags.getBirthCells();
    const auto& death_cells = tags.getDeathCells();
    const auto& birth_times = tags.getBirthTimes();
    const auto& death_times = tags.getDeathTimes();
    
    // Store generating sets for conflict computation
    std::vector<int> real_birth_cells;
    real_birth_cells.reserve(birth_cells.size());
    std::vector<int> real_death_cells;
    real_death_cells.reserve(death_cells.size());
    std::vector<std::vector<int>> real_generating_sets_fore_indices;
    real_generating_sets_fore_indices.reserve(birth_cells.size());
    std::vector<std::vector<int>> real_generating_sets_back_indices;
    real_generating_sets_back_indices.reserve(death_cells.size());

    
    
    std::cout << "\nGenerating sets for birth cells:" << std::endl;
    
    // Create vector masks for faster lookup (O(1) direct access instead of hash lookup)
    int max_index = std::max(boundary_matrix.rows(), boundary_matrix.cols());
    std::vector<bool> remaining_negative_mask(max_index, false);
    std::vector<bool> remaining_positive_mask(max_index, false);
    
    for (int idx : remaining_negative) {
        if (idx >= 0 && idx < max_index) {
            remaining_negative_mask[idx] = true;
        }
    }
    for (int idx : remaining_positive) {
        if (idx >= 0 && idx < max_index) {
            remaining_positive_mask[idx] = true;
        }
    }

    // Skip first topK pairs when computing generating sets
    size_t start_idx = std::min(static_cast<size_t>(topK), birth_cells.size());
    std::cout << "Using " << (birth_cells.size() - start_idx) << " pairs after ignoring top " << topK << " pairs\n\n";
    
    // Alternative version: find first pair where birth_time + death_time < threshold
    /*
    size_t start_idx = 0;
    for (size_t i = 0; i < birth_cells.size(); ++i) {
        if (birth_times[i] + death_times[i] < 0.95) {
            start_idx = i;
            break;
        }
    }
    std::cout << "Using " << (birth_cells.size() - start_idx) << " pairs after finding first pair with birth_time + death_time < 0.8 (start_idx = " << start_idx << ")\n\n";
    */
    std::vector<std::vector<int>> fore_col_lists(tags.getForegroundMatrix().cols());
    
    // Pre-allocate space in lists
    for (auto& list : fore_col_lists) list.reserve(6);
    
    // Build lists in a single pass
    for (int k = 0; k < tags.getForegroundMatrix().outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(tags.getForegroundMatrix(), k); it; ++it) {
            int row = it.row();
            int col = it.col();
            fore_col_lists[row].push_back(col);
        }
    }

    std::vector<std::vector<int>> back_row_lists(tags.getBackgroundMatrix().rows());
    
    // Pre-allocate space in lists
    for (auto& list : back_row_lists) list.reserve(6);
    
    int count_core = 0;
    int count_neighborhood = 0;
    // Build lists in a single pass
    for (int k = 0; k < tags.getBackgroundMatrix().outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(tags.getBackgroundMatrix(), k); it; ++it) {
            int row = it.row();
            int col = it.col();
            back_row_lists[col].push_back(row);
        }
    }
    
    auto dimensions = boundary_matrix.getDimensions();

    std::vector<int> old_birth_indices;
    old_birth_indices.reserve(birth_cells.size());
    std::vector<int> old_death_indices;
    old_death_indices.reserve(death_cells.size());
    std::vector<int> new_birth_indices;
    new_birth_indices.reserve(death_cells.size());

    
    //TESTs
    int cavity_counter = 0;
    int handle_counter = 0;
    int component_counter = 0;
    /*
    cavity_counter = 0; 

    for (size_t i = 0; i < birth_cells.size(); ++i) {
        if (cavity_counter < 3 && dimensions[birth_cells[i]] == 2) {
            cavity_counter++;
            
            auto curr_generator = sorted_pairs[i].generator;
            std::cout << "CAVITY birth time: " << birth_times[i] << std::endl;
            std::cout << "CAVITY death time: " << death_times[i] << std::endl;

            //output the generators into a file
            std::string generators_file_path = (topomc_output_root / "generators" / (std::to_string(cavity_counter) + ".txt")).string();
            std::ofstream generators_file(generators_file_path);
            if (!generators_file.is_open()) {
                std::cerr << "Failed to open generators file" << std::endl;
                return 1;
            }
            for (int generator : curr_generator) {
                generators_file << generator << std::endl;
            
            continue;
            }
        }
    }
    */
    std::vector<bool> isolated_fore_vec(boundary_matrix.rows(), true);

    for (size_t i = 0; i< remaining_negative.size(); ++i){
        int idx = remaining_negative[i];
        std::vector<int> col_indices;
        boundary_matrix.getColIndices(idx, col_indices);

        for (int idx2 : col_indices){
            isolated_fore_vec[idx2] = false;
        }
    }

    std::vector<std::vector<int>> all_generating_sets_fore_indices_test;
    all_generating_sets_fore_indices_test.reserve(birth_cells.size());
    std::vector<std::vector<int>> all_generating_sets_back_indices_test;
    all_generating_sets_back_indices_test.reserve(death_cells.size());
    
    for (size_t i = start_idx; i < birth_cells.size(); ++i) {

        if (cavity_counter < cavitySkip && dimensions[birth_cells[i]] == 2){
            //std::cout << "Skipping feature: " << birth_cells[i] << " (alpha=" << birth_times[i] << ")" << std::endl;
            //std::cout << "Death cell: " << death_cells[i] << " (alpha=" << death_times[i] << ")" << std::endl;
            cavity_counter++;
            continue;
        }

        if (handle_counter < handleSkip && dimensions[birth_cells[i]] == 1){
            //std::cout << "Skipping feature: " << birth_cells[i] << " (alpha=" << birth_times[i] << ")" << std::endl;
            //std::cout << "Death cell: " << death_cells[i] << " (alpha=" << death_times[i] << ")" << std::endl;
            handle_counter++;
            continue;
        }

        if (component_counter < componentSkip && dimensions[birth_cells[i]] == 0){
            //std::cout << "Skipping feature: " << birth_cells[i] << " (alpha=" << birth_times[i] << ")" << std::endl;
            //std::cout << "Death cell: " << death_cells[i] << " (alpha=" << death_times[i] << ")" << std::endl;
            component_counter++;
            continue;
        }


        int birth_cell = birth_cells[i];
        int death_cell = death_cells[i];
        double birth_time = birth_times[i];
        double death_time = death_times[i];

        bool isolated_fore = isolated_fore_vec[birth_cell];

        std::vector<int> fore_set;
        if (isolated_fore) {
            //fore_set = tags.getGeneratingSetForeFast(birth_cell);
            fore_set = tags.getGeneratingSetForeFast2(fore_col_lists, birth_cell);
        } // else leave as empty set

        // Check isolation for BackFast - using pre-computed full_row_lists
        bool isolated_back = true;
        std::vector<int> col_indices;
        boundary_matrix.getColIndices(death_cell, col_indices);
        for (int idx : col_indices) {
            if (idx >= 0 && idx < static_cast<int>(remaining_positive_mask.size()) && remaining_positive_mask[idx]) {
                isolated_back = false;
                break;
            }
        }

        std::vector<int> back_set;
        if (isolated_back) {
            //back_set = tags.getGeneratingSetBackFast(death_cell);
            back_set = tags.getGeneratingSetBackFast2(back_row_lists, death_cell);
        } // else leave as empty set

        // Store foreground if non-empty
        if (!fore_set.empty()&& birth_time>=core) {
            real_birth_cells.push_back(birth_cell);
            real_generating_sets_fore_indices.push_back(fore_set);
            old_birth_indices.push_back(i);
        }
        all_generating_sets_fore_indices_test.push_back(fore_set);
        //count the number of birth_times>=core
        if (birth_time<core && !fore_set.empty()) {
            count_core++;
        }

        if (!back_set.empty()&& death_time<=neighborhood) {
            real_death_cells.push_back(death_cell);
            real_generating_sets_back_indices.push_back(back_set);
            old_death_indices.push_back(i);
        }
        all_generating_sets_back_indices_test.push_back(back_set);

        if (death_time>neighborhood && !back_set.empty()) {
            count_neighborhood++;
        }

    }
    std::cout << "(invalid birth) Number of birth_times<core and fore_set is not empty: " << count_core << std::endl;
    std::cout << "(invalid death)Number of death_times>neighborhood and back_set is not empty: " << count_neighborhood << std::endl;

    // Free memory from the lists

    fore_col_lists.clear();
    fore_col_lists.shrink_to_fit();
    back_row_lists.clear();
    back_row_lists.shrink_to_fit();

    /*
    
    std::vector<int> real_birth_cells_test;
    real_birth_cells_test.reserve(1000);
    std::vector<int> real_death_cells_test;
    real_death_cells_test.reserve(1000);
    std::vector<std::vector<int>> real_generating_sets_fore_indices_test;
    real_generating_sets_fore_indices_test.reserve(1000);
    std::vector<std::vector<int>> real_generating_sets_back_indices_test;
    real_generating_sets_back_indices_test.reserve(1000);

    std::cout << "Starting test for i,j,k of 0-10: " << birth_cells.size() << std::endl;
    for (int i=0; i<1; i++){
        std::cout << "i=" << i << std::endl;
        for (int j=0; j<456; j++){
            for (int k=0; k<1; k++){
                real_birth_cells_test.clear();
                real_death_cells_test.clear();
                real_generating_sets_fore_indices_test.clear();
                real_generating_sets_back_indices_test.clear();

                real_birth_cells_test.reserve(1000);
                real_death_cells_test.reserve(1000);
                real_generating_sets_fore_indices_test.reserve(1000);
                real_generating_sets_back_indices_test.reserve(1000);

                //ignore top i,j,k in its corresponding dimension
                int iCount = 0;
                int jCount = 0;
                int kCount = 0;

                //populate test generating sets
                for (size_t l=0; l<birth_cells.size(); l++){

                    if (dimensions[birth_cells[l]] == 0 && iCount < i){
                        iCount++;
                        continue;
                    }
                    if (dimensions[birth_cells[l]] == 1 && jCount < j){
                        jCount++;
                        continue;
                    }
                    if (dimensions[birth_cells[l]] == 2 && kCount < k){
                        kCount++;
                        continue;
                    }
                    
                    if(!all_generating_sets_fore_indices_test[l].empty()){
                        real_birth_cells_test.push_back(birth_cells[l]);
                        real_generating_sets_fore_indices_test.push_back(all_generating_sets_fore_indices_test[l]);
                    }
                    if(!all_generating_sets_back_indices_test[l].empty()){
                        real_death_cells_test.push_back(death_cells[l]);
                        real_generating_sets_back_indices_test.push_back(all_generating_sets_back_indices_test[l]);
                    }
                }


                //compute conflicts
                auto conflicts_test = tags.computeConflicts(real_birth_cells_test, real_death_cells_test,
                                         real_generating_sets_fore_indices_test,
                                         real_generating_sets_back_indices_test);
                
                auto [cut_birth_nodes_test, cut_death_nodes_test] = compute_min_cut(
                    real_birth_cells_test,
                    real_death_cells_test,
                    real_generating_sets_fore_indices_test,
                    real_generating_sets_back_indices_test,
                    conflicts_test,
                    alphas
                );

                                // Compute ISET_birth and ISET_death (complements)
                std::vector<int> ISET_birth_test = compute_iset_indices(cut_birth_nodes_test, static_cast<int>(real_birth_cells_test.size()));
                std::vector<int> ISET_death_test = compute_iset_indices(cut_death_nodes_test, static_cast<int>(real_death_cells_test.size()));

                int remaining_topology = birth_cells.size() -i-j-k - ISET_birth_test.size() - ISET_death_test.size();
                if (remaining_topology!= 0){

                    std::cout << "MinCut didn't remove all topology for i=" << i << ", j=" << j << ", k=" << k << std::endl;
                    std::cout << "Remaining topology: " << remaining_topology << std::endl;

                }
  
            }
        }
    }

    if (real_birth_cells_test.size() < 10000){
        std::cout << "End of testing!, returning 0: " << real_birth_cells_test.size() << std::endl;
        return 0;
    }
    */
    // Compute conflicts
    auto conflicts = tags.computeConflicts(real_birth_cells, real_death_cells,
                                         real_generating_sets_fore_indices,
                                         real_generating_sets_back_indices);

    // Compute min-cut using the new function
    
    auto [cut_birth_nodes, cut_death_nodes] = compute_min_cut(
        real_birth_cells,
        real_death_cells,
        real_generating_sets_fore_indices,
        real_generating_sets_back_indices,
        conflicts,
        alphas
    );

    // Compute ISET_birth and ISET_death (complements)
    std::vector<int> ISET_birth = compute_iset_indices(cut_birth_nodes, static_cast<int>(real_birth_cells.size()));
    std::vector<int> ISET_death = compute_iset_indices(cut_death_nodes, static_cast<int>(real_death_cells.size()));

    std::vector<int> cut_fore_cells_vec = union_generating_sets(real_generating_sets_fore_indices, ISET_birth);
    std::vector<int> fill_back_cells_vec = union_generating_sets(real_generating_sets_back_indices, ISET_death);

    
    std::cout << "Number of cut_fore_cells_vec: " << ISET_birth.size() << std::endl;
    std::cout << "Number of fill_back_cells_vec: " << ISET_death.size() << std::endl;
    std::cout << "Number of total used cells: " << ISET_birth.size() + ISET_death.size() << std::endl;

    
    //lets actually find the used indices (for either birth or death) and eventually find the complement (unused indices)
    std::vector<int> used_indices;
    used_indices.reserve(old_birth_indices.size() + old_death_indices.size());
    
    // Collect all old indices that have either a birth or death used
    for (int idx : ISET_birth) {
        if (idx < old_birth_indices.size()) {
            used_indices.push_back(old_birth_indices[idx]);
        }
    }
    for (int idx : ISET_death) {
        if (idx < old_death_indices.size()) {
            used_indices.push_back(old_death_indices[idx]);
        }
    }
    
    // Remove duplicates and sort
    std::sort(used_indices.begin(), used_indices.end());
    used_indices.erase(std::unique(used_indices.begin(), used_indices.end()), used_indices.end());
    
    // Find complement (unused indices)
    std::vector<int> unused_indices;
    int total_indices = static_cast<int>(birth_cells.size());
    for (int i = 0; i < total_indices; ++i) {
        if (std::find(used_indices.begin(), used_indices.end(), i) == used_indices.end()) {
            unused_indices.push_back(i);
        }
    }
    
    std::cout << "Number of used indices: " << used_indices.size() << std::endl;
    std::cout << "Number of unused indices: " << unused_indices.size() << std::endl;

    //create a list of certain birth and death cells that are used for the solution
    std::vector<int> used_birth_cells;
    used_birth_cells.reserve(used_indices.size());
    std::vector<int> used_death_cells;
    used_death_cells.reserve(used_indices.size());
    for (int idx : used_indices) {
        used_birth_cells.push_back(birth_cells[idx]);
        used_death_cells.push_back(death_cells[idx]);
    }
    
    


    // === FINAL TASK: Update alphas, reorder, and build new boundary matrix ===
    
    // 1. Set new alpha values for cut/filled cells
    std::vector<double> new_alphas = alphas;

    //find the highest negative alpha
    double highest_negative_alpha = -std::numeric_limits<double>::infinity();
    for (double alpha : alphas) {
        if (alpha < 0.0 && alpha > highest_negative_alpha) {
            highest_negative_alpha = alpha;
        }
    }
    std::cout << "Highest negative alpha: " << highest_negative_alpha << std::endl;

    //find the lowest positive alpha
    double lowest_positive_alpha = std::numeric_limits<double>::infinity();
    for (double alpha : alphas) {
        if (alpha > 0.0 && alpha < lowest_positive_alpha) {
            lowest_positive_alpha = alpha;
        }
    }
    std::cout << "Lowest positive alpha: " << lowest_positive_alpha << std::endl;
    
    for (int idx : cut_fore_cells_vec) {
        new_alphas[idx] = lowest_positive_alpha/2.0;
        new_alphas[idx] = 0.001;
    }
    for (int idx : fill_back_cells_vec) {
        new_alphas[idx] = highest_negative_alpha/2.0;
        new_alphas[idx] = -0.001;
    }
    
    // Output new alphas to global output directory first
    /*
    std::ofstream global_alpha_file(topomc_output_root / "alphas_updated.txt");
    if (!global_alpha_file.is_open()) {
        std::cerr << "Failed to open global alpha output file" << std::endl;
        return 1;
    }
    //global_alpha_file << new_alphas.size() << std::endl;
    for (double alpha : new_alphas) {
        if (std::isinf(alpha)) {
            global_alpha_file << "infinity" << std::endl;
        } else {
            global_alpha_file << alpha << std::endl;
        }
    }
    global_alpha_file.close();
    */

    // 2. Compute stable sort permutation (new_to_old)
    std::vector<int> new_to_old(new_alphas.size());
    std::iota(new_to_old.begin(), new_to_old.end(), 0);
    std::stable_sort(new_to_old.begin(), new_to_old.end(),
        [&](int i, int j) {
            if (new_alphas[i] != new_alphas[j])
                return new_alphas[i] < new_alphas[j];
            return i < j;
        }
    );

    // 3. Build old_to_new index map
    std::vector<int> old_to_new(new_alphas.size());
    for (size_t i = 0; i < new_to_old.size(); ++i) {
        old_to_new[new_to_old[i]] = i;
    }

    // 4. Build new boundary matrix
    Eigen::SparseMatrix<int> new_matrix(new_alphas.size(), new_alphas.size());
    
    // Pre-allocate triplets with exact size from original matrix
    std::vector<Eigen::Triplet<int>> triplets;
    triplets.reserve(boundary_matrix.getMatrix().nonZeros());
    
    // Process all columns at once
    for (int new_col = 0; new_col < (int)new_alphas.size(); ++new_col) {
        int old_col = new_to_old[new_col];
        std::vector<int> old_rows;
        boundary_matrix.getColIndices(old_col, old_rows);
        for (int old_row : old_rows) {
            triplets.emplace_back(old_to_new[old_row], new_col, 1);
        }
    }
    new_matrix.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << std::endl;
    std::cout << "New boundary matrix nonzeros: " << new_matrix.nonZeros() << std::endl;

    // 5. Update final alphas and dimensions according to new order
    std::vector<double> final_alphas(new_alphas.size());
    std::vector<int> final_dimensions(new_alphas.size());
    const auto& old_dimensions = boundary_matrix.getDimensions();
    
    for (size_t i = 0; i < new_alphas.size(); ++i) {
        final_alphas[i] = new_alphas[new_to_old[i]];
        final_dimensions[i] = old_dimensions[new_to_old[i]];
    }

    
    // Output final alphas to file
    /*
    {
        std::ofstream alpha_file(topomc_output_root / "new_alphas.txt");
        if (!alpha_file.is_open()) {
            std::cerr << "Failed to open alpha output file" << std::endl;
            return 1;
        }
        alpha_file << final_alphas.size() << std::endl;
        for (double alpha : final_alphas) {
            if (std::isinf(alpha)) {
                alpha_file << "infinity" << std::endl;
            } else {
                alpha_file << alpha << std::endl;
            }
        }
        alpha_file.close();
    } // alpha_file is automatically closed here
    */
    /*
    // Output birth cells
    {
        std::ofstream birth_file(topomc_output_root / "birth_cells.txt");
        if (!birth_file.is_open()) {
            std::cerr << "Failed to open birth cells output file" << std::endl;
            return 1;
        }
        for (int cell : real_birth_cells) {
            birth_file << cell << std::endl;
        }
        birth_file.close();
    }

    // Output death cells
    {
        std::ofstream death_file(topomc_output_root / "death_cells.txt");
        if (!death_file.is_open()) {
            std::cerr << "Failed to open death cells output file" << std::endl;
            return 1;
        }
        for (int cell : real_death_cells) {
            death_file << cell << std::endl;
        }
        death_file.close();
    }
    */
    // Output fore generating sets
    /*
    {
        std::ofstream fore_sets_file(topomc_output_root / "fore_generating_sets.txt");
        if (!fore_sets_file.is_open()) {
            std::cerr << "Failed to open fore generating sets output file" << std::endl;
            return 1;
        }
        for (const auto& generator_set : real_generating_sets_fore_indices) {
            for (int cell : generator_set) {
                fore_sets_file << cell << std::endl;
            }
        }
        fore_sets_file.close();
    }
        */
    /*
    // Output fore generators
    {
        std::ofstream fore_file(topomc_output_root / "fore_generators.txt");
        if (!fore_file.is_open()) {
            std::cerr << "Failed to open fore generators output file" << std::endl;
            return 1;
        }
        for (const auto& pair : sorted_pairs) {
            for (int cell : pair.generator) {
                fore_file << cell << std::endl;
            }
        }
        fore_file.close();
        
    }

    // Output fore values (using absolute birth values)
    {
        std::ofstream fore_values_file(topomc_output_root / "fore_values.txt");
        if (!fore_values_file.is_open()) {
            std::cerr << "Failed to open fore values output file" << std::endl;
            return 1;
        }
        // For each generating set
        for (size_t i = 0; i < real_generating_sets_fore_indices.size(); ++i) {
            int birth_cell = real_birth_cells[i];
            double birth_value = std::abs(alphas[birth_cell]);
            // Write the same birth value for each index in this generating set
            const auto& generator_set = real_generating_sets_fore_indices[i];
            for (size_t j = 0; j < generator_set.size(); ++j) {
                fore_values_file << birth_value << std::endl;
            }
        }
        fore_values_file.close();
    }
    */
    // Output back generating sets

    vector<int> removal_order;
    removal_order.reserve(alphas.size());
    for (size_t i = 0; i < alphas.size(); ++i) {
        removal_order.push_back(0);
    }

    for (size_t i=0; i< background_removed_pairs.size(); ++i) {
        removal_order[background_removed_pairs[i].second] = i;
        removal_order[background_removed_pairs[i].first] = i;
    
    }
    /*
    {
        std::ofstream back_sets_file(topomc_output_root / "back_generating_sets.txt");
        std::ofstream back_sets_file_order(topomc_output_root / "back_generating_sets_order.txt");
        if (!back_sets_file.is_open()) {
            std::cerr << "Failed to open back generating sets output file" << std::endl;
            return 1;
        }
        if (!back_sets_file_order.is_open()) {
            std::cerr << "Failed to open back generating sets order output file" << std::endl;
            return 1;
        }
        for (const auto& generator_set : real_generating_sets_back_indices) {
            for (int cell : generator_set) {
                back_sets_file << cell << std::endl;
                back_sets_file_order << removal_order[cell] << std::endl;
            }
        }
        back_sets_file.close();
        back_sets_file_order.close();
    }
*/
    /*
    // Output back values (using death cell values)
    {
        std::ofstream back_values_file(topomc_output_root / "back_values.txt");
        if (!back_values_file.is_open()) {
            std::cerr << "Failed to open back values output file" << std::endl;
            return 1;
        }
        // For each generating set
        for (size_t i = 0; i < real_generating_sets_back_indices.size(); ++i) {
            int death_cell = real_death_cells[i];
            double death_value = std::abs(alphas[death_cell]);
            // Write the same death value for each index in this generating set
            const auto& generator_set = real_generating_sets_back_indices[i];
            for (size_t j = 0; j < generator_set.size(); ++j) {
                back_values_file << death_value << std::endl;
            }
        }
        back_values_file.close();
    }

     */
     /*
    // Output dual generators
    {
        std::ofstream dual_file(topomc_output_root / "dual_generators.txt");
        if (!dual_file.is_open()) {
            std::cerr << "Failed to open dual generators output file" << std::endl;
            return 1;
        }
        // Output all entries except the last one
        for (size_t i = 0; i < protected_indices_background_vec.size() - 1; ++i) {
            dual_file << protected_indices_background_vec[i] << std::endl;
        }
        dual_file.close();
    }

    // Output ISET birth generating sets
    {
        std::ofstream iset_birth_file(topomc_output_root / "iset_birth_generating_sets.txt");
        if (!iset_birth_file.is_open()) {
            std::cerr << "Failed to open ISET birth generating sets output file" << std::endl;
            return 1;
        }
        for (int birth_idx : ISET_birth) {
            // Get the corresponding generating set
            const auto& generator_set = real_generating_sets_fore_indices[birth_idx];
            for (int cell : generator_set) {
                iset_birth_file << cell << std::endl;
            }
        }
        iset_birth_file.close();
    }

    // Output ISET death generating sets
    {
        std::ofstream iset_death_file(topomc_output_root / "iset_death_generating_sets.txt");
        if (!iset_death_file.is_open()) {
            std::cerr << "Failed to open ISET death generating sets output file" << std::endl;
            return 1;
        }
        for (int death_idx : ISET_death) {
            // Get the corresponding generating set
            const auto& generator_set = real_generating_sets_back_indices[death_idx];
            for (int cell : generator_set) {
                iset_death_file << cell << std::endl;
            }
        }
        iset_death_file.close();
    }
    */

    // Stop when topology is cleared or iteration budget is exhausted; otherwise reload in memory for next round.
    long long remaining_topology =
        static_cast<long long>(birth_cells.size()) -
        static_cast<long long>(cavity_counter + handle_counter + component_counter) -
        static_cast<long long>(ISET_birth.size() + ISET_death.size());

    const bool done_topo = (remaining_topology <= 0);
    const bool hit_cap = (iteration + 1 >= max_topo_iters);

    std::cout << "Remaining topology count (after cut): " << remaining_topology << std::endl;

    if (outs) {
        outs->remaining_negative.clear();
        outs->protected_indices.clear();
        for (int idx : remaining_negative) {
            if (alphas[idx] >= core) {
                const int o = col_original_id[static_cast<size_t>(idx)];
                outs->remaining_negative.push_back(o);
                outs->protected_indices.push_back(o);
            }
        }
    }

    if (done_topo || hit_cap) {
        if (hit_cap && !done_topo) {
            std::cerr << "Warning: stopping TopoMinCut after maxTopologyIterations="
                      << max_topo_iters << " with remaining_topology=" << remaining_topology << std::endl;
        }
        if (outs) {
            outs->alphas_updated.assign(static_cast<size_t>(n0), 0.0);
            if (done_topo && iteration == 0) {
                for (size_t pos = 0; pos < new_alphas.size(); ++pos) {
                    outs->alphas_updated[static_cast<size_t>(col_original_id[pos])] = new_alphas[pos];
                }
            } else {
                std::vector<int> col_after(new_to_old.size());
                for (size_t nc = 0; nc < new_to_old.size(); ++nc) {
                    col_after[nc] = col_original_id[static_cast<size_t>(new_to_old[nc])];
                }
                for (size_t nc = 0; nc < final_alphas.size(); ++nc) {
                    outs->alphas_updated[static_cast<size_t>(col_after[nc])] = final_alphas[nc];
                }
            }
            outs->permuted_boundary_matrix = new_matrix;
            outs->completed_second_iteration = (iteration >= 1);
            outs->iterations_completed = iteration + 1;
        }
        return 0;
    }

    std::cout << "Topology remains; preparing iteration " << (iteration + 2) << " (in-memory reload)." << std::endl;

    {
        std::vector<int> col_after(new_to_old.size());
        for (size_t nc = 0; nc < new_to_old.size(); ++nc) {
            col_after[nc] = col_original_id[static_cast<size_t>(new_to_old[nc])];
        }
        col_original_id.swap(col_after);
    }

    if (!boundary_matrix.loadFromEigenCopy(std::move(new_matrix), final_dimensions)) {
        std::cerr << "Failed to reload boundary matrix from Eigen for next iteration" << std::endl;
        return 1;
    }
    if (!boundary_matrix.loadAlphasVector(final_alphas)) {
        std::cerr << "Failed to reload alphas for next iteration" << std::endl;
        return 1;
    }
    }

}

int runFromEigenSparse(const Eigen::SparseMatrix<int>& matrix,
                       const std::vector<int>& dimensions,
                       const std::vector<double>& alphas,
                       const RunParams& params,
                       RunOutputs* outs) {
    BoundaryMatrix bm;
    Eigen::SparseMatrix<int> copy = matrix;
    if (!bm.loadFromEigenCopy(std::move(copy), dimensions)) {
        return 1;
    }
    if (!bm.loadAlphasVector(alphas)) {
        return 1;
    }
    return runFromBoundaryMatrix(bm, params, outs);
}

} // namespace topomincut 