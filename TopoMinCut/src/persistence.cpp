#include "persistence.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

Persistence::Persistence(const BoundaryMatrix& boundary_matrix, const std::string& matrix_file)
    : boundary_matrix_(boundary_matrix), matrix_file_(matrix_file) {
    // Initialize matrices with same dimensions as boundary matrix
    reduced_matrix_ = Eigen::SparseMatrix<int>(boundary_matrix.rows(), boundary_matrix.cols());
    dual_reduced_matrix_ = Eigen::SparseMatrix<int>(boundary_matrix.rows(), boundary_matrix.cols());
}

Persistence::~Persistence() {}

void Persistence::computePersistence(double core, double neighborhood) {
    try {
        std::cout << "\n=== Starting Persistence Computation ===" << std::endl;
        
        // Step 1: Load PHAT boundary matrix from in-memory Eigen boundary matrix (PHAT matches file layout).
        std::cout << "\nStep 1: Loading matrix into PHAT from BoundaryMatrix..." << std::endl;
        const bool verbose_compare = (boundary_matrix_.cols() <= 500);
        loadFromBoundaryMatrix(verbose_compare);
        std::cout << "Matrix loaded successfully. Dimensions: " << phat_matrix_.get_num_cols() << std::endl;
        
        // Step 2: Perform exhaustive reduction
        std::cout << "\nStep 2: Starting exhaustive reduction..." << std::endl;
        performExhaustiveReduction();
        std::cout << "Exhaustive reduction completed. Found " << birth_cells_.size() << " pairs." << std::endl;
        
        // Step 3: Convert to Eigen format (skip for huge matrices to save memory)
        std::cout << "\nStep 3: Converting to Eigen format..." << std::endl;
        if (boundary_matrix_.cols() <= 100) {
            convertFromPHAT();
            std::cout << "Conversion completed. Matrix nonzeros: " << reduced_matrix_.nonZeros() << std::endl;
        } else {
            std::cout << "Skipping Eigen conversion for large matrix (" << boundary_matrix_.cols() << ")" << std::endl;
        }

        // Step 4: Filter pairs based on alpha values
        std::cout << "\nStep 4: Filtering pairs based on alpha values..." << std::endl;
        filterPairs(core, neighborhood);
        std::cout << "Filtering completed. Found " << filtered_pairs_.size() << " filtered pairs." << std::endl;
        
        // Step 5: Perform dual reduction
        std::cout << "\nStep 5: Starting dual reduction..." << std::endl;
        performDualReduction();
        std::cout << "Dual reduction completed." << std::endl;

        
        std::cout << "\n=== Persistence Computation Completed ===" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error in computePersistence: " << e.what() << std::endl;
        throw;
    }
}

void Persistence::loadFromBoundaryMatrix(bool verbose_compare_with_eigen) {
    const int number_of_columns = boundary_matrix_.cols();
    if (number_of_columns <= 0) {
        throw std::runtime_error("Boundary matrix has no columns");
    }

    phat_matrix_.set_num_cols(number_of_columns);
    dual_phat_matrix_.set_num_cols(number_of_columns);

    const auto& dims = boundary_matrix_.getDimensions();
    const Eigen::SparseMatrix<int>& eigen_mat = boundary_matrix_.getMatrix();

    for (int col = 0; col < number_of_columns; ++col) {
        phat_matrix_.set_dim(col, static_cast<phat::dimension>(dims[col]));
        dual_phat_matrix_.set_dim(col, static_cast<phat::dimension>(dims[col]));

        phat::column temp_col;
        for (Eigen::SparseMatrix<int>::InnerIterator it(eigen_mat, col); it; ++it) {
            if (it.value() != 0) {
                temp_col.push_back(static_cast<phat::index>(it.row()));
            }
        }
        std::sort(temp_col.begin(), temp_col.end());
        phat_matrix_.set_col(col, temp_col);
        dual_phat_matrix_.set_col(col, temp_col);
    }

    if (!verbose_compare_with_eigen) {
        std::cout << "(Skipping Eigen vs PHAT column verify for large matrix.)\n";
        return;
    }

    int mismatched_columns = 0;
    for (int c = 0; c < number_of_columns; c++) {
        phat::column phat_col;
        phat_matrix_.get_col(c, phat_col);
        std::sort(phat_col.begin(), phat_col.end());

        std::vector<int> eigen_col;
        boundary_matrix_.getColIndices(c, eigen_col);
        std::sort(eigen_col.begin(), eigen_col.end());

        if (phat_col.size() != eigen_col.size()) {
            mismatched_columns++;
            continue;
        }
        bool mismatch = false;
        for (size_t i = 0; i < phat_col.size(); i++) {
            if (phat_col[i] != eigen_col[i]) {
                mismatch = true;
                break;
            }
        }
        if (mismatch) {
            mismatched_columns++;
        }
    }
    std::cout << "Eigen vs PHAT column mismatch count: " << mismatched_columns << " / " << number_of_columns << std::endl;
}

void Persistence::loadFromFile(const std::string& filename) {
    // First count number of columns like PHAT does
    std::string cur_line;
    std::ifstream dummy(filename.c_str());
    if (dummy.fail()) {
        throw std::runtime_error("Could not open matrix file: " + filename);
    }

    // Skip first line (dimension info)
    getline(dummy, cur_line);

    int number_of_columns = 0;
    while (getline(dummy, cur_line)) {
        cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (cur_line != "" && cur_line[0] != '#') {
            number_of_columns++;
        }
    }
    dummy.close();

    // Set dimensions in PHAT matrices
    phat_matrix_.set_num_cols(number_of_columns);
    dual_phat_matrix_.set_num_cols(number_of_columns);

    // Read matrix entries like PHAT does
    std::ifstream input_stream(filename.c_str());
    if (input_stream.fail()) {
        throw std::runtime_error("Could not open matrix file for reading: " + filename);
    }

    // Skip first line (dimension info)
    getline(input_stream, cur_line);

    phat::column temp_col;
    int cur_col = -1;
    while (getline(input_stream, cur_line)) {
        cur_line.erase(cur_line.find_last_not_of(" \t\n\r\f\v") + 1);
        if (cur_line != "" && cur_line[0] != '#') {
            cur_col++;
            std::stringstream ss(cur_line);
            
            // Read dimension (first number in line)
            int64_t temp_dim;
            ss >> temp_dim;
            phat_matrix_.set_dim(cur_col, (phat::dimension)temp_dim);
            dual_phat_matrix_.set_dim(cur_col, (phat::dimension)temp_dim);

            // Read column entries
            temp_col.clear();
            int64_t temp_index;
            while (ss.good()) {
                ss >> temp_index;
                temp_col.push_back((phat::index)temp_index);
            }
            
            // Sort column entries like PHAT does
            std::sort(temp_col.begin(), temp_col.end());
            phat_matrix_.set_col(cur_col, temp_col);
            dual_phat_matrix_.set_col(cur_col, temp_col);
        }
    }
    input_stream.close();

    // Compare non-zero entries between Eigen and PHAT matrices
    std::cout << "\nComparing non-zero entries between Eigen and PHAT matrices:" << std::endl;
    int mismatched_columns = 0;
    for (int c = 0; c < number_of_columns; c++) {
        // Get PHAT column entries
        phat::column phat_col;
        phat_matrix_.get_col(c, phat_col);
        std::sort(phat_col.begin(), phat_col.end());

        // Get Eigen column entries
        std::vector<int> eigen_col;
        boundary_matrix_.getColIndices(c, eigen_col);
        std::sort(eigen_col.begin(), eigen_col.end());

        // Compare
        if (phat_col.size() != eigen_col.size()) {
            mismatched_columns++;
            continue;
        }

        // Check if entries match
        bool mismatch = false;
        for (size_t i = 0; i < phat_col.size(); i++) {
            if (phat_col[i] != eigen_col[i]) {
                mismatch = true;
                break;
            }
        }
        if (mismatch) {
            mismatched_columns++;
        }
    }

    std::cout << "Total columns with mismatches: " << mismatched_columns << " out of " << number_of_columns << std::endl;

    // Print total non-zeros for both matrices
    std::cout << "\nMatrix Statistics:" << std::endl;
    
    // 1. Eigen matrix non-zeros (from boundary_matrix_)
    std::cout << "Eigen Matrix (from boundary_matrix_):" << std::endl;
    std::cout << "Total non-zeros: " << boundary_matrix_.getMatrix().nonZeros() << std::endl;

    // 2. PHAT matrix non-zeros
    std::cout << "\nPHAT Matrix:" << std::endl;
    int total_non_zeros = 0;
    for (int col = 0; col < number_of_columns; ++col) {
        phat::column col_data;
        phat_matrix_.get_col(col, col_data);
        total_non_zeros += col_data.size();
    }
    std::cout << "Total non-zeros: " << total_non_zeros << std::endl;
}

void Persistence::performExhaustiveReduction() {
    try {
        std::cout << "  Starting PHAT reduction..." << std::endl;
        phat::persistence_pairs pairs;
        
        std::cout << "  Computing persistence pairs..." << std::endl;
        phat::compute_persistence_pairs<phat::exhaustive_compress_reduction>(pairs, phat_matrix_);
        

        //TESTING STANDARD REDUCTION
        //phat::compute_persistence_pairs<phat::standard_reduction>(pairs, phat_matrix_);

        // Pre-allocate vectors with known size
        size_t num_pairs = pairs.get_num_pairs();
        std::cout << "  Found " << num_pairs << " total pairs" << std::endl;
        
        birth_cells_.reserve(num_pairs);
        death_cells_.reserve(num_pairs);
        
        // Clear existing data
        birth_cells_.clear();
        death_cells_.clear();
        
        // Process pairs
        std::cout << "  Processing pairs..." << std::endl;
        size_t step_size = (num_pairs > 10) ? (num_pairs / 10) : 1;  // Show ~10 progress updates
        for (int i = 0; i < num_pairs; ++i) {
            if (i % step_size == 0) {
                std::cout << "    Processed " << i << " pairs..." << std::endl;
            }
            auto pair = pairs.get_pair(i);
            birth_cells_.push_back(pair.first);
            death_cells_.push_back(pair.second);
        }
        std::cout << "  All pairs processed." << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in performExhaustiveReduction: " << e.what() << std::endl;
        throw;
    }
}

void Persistence::performDualReduction() {
    try {
        std::cout << "  Starting dual reduction..." << std::endl;
        
        // Dualize the matrix
        std::cout << "  Dualizing matrix..." << std::endl;
        phat::dualize(dual_phat_matrix_);
        
        // Perform exhaustive reduction on dualized matrix
        std::cout << "  Computing dual persistence pairs..." << std::endl;
        phat::persistence_pairs dual_pairs;

        phat::compute_persistence_pairs<phat::exhaustive_compress_reduction>(dual_pairs, dual_phat_matrix_);
        //TESTING
        //phat::compute_persistence_pairs<phat::standard_reduction>(dual_pairs, dual_phat_matrix_);

        std::cout << "  Dual reduction completed. Found " << dual_pairs.get_num_pairs() << " dual pairs" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error in performDualReduction: " << e.what() << std::endl;
        throw;
    }
}

void Persistence::filterPairs(double core, double neighborhood) {
    std::cout << "  Starting pair filtering..." << std::endl;
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    std::cout << "  Got " << alphas.size() << " alpha values" << std::endl;

    filtered_pairs_.clear();
    filtered_pairs_.reserve(birth_cells_.size());
    std::cout << "  Processing " << birth_cells_.size() << " pairs for filtering" << std::endl;

    auto cell_dimensions = boundary_matrix_.getDimensions();
    int component_counter = 0;
    int handle_counter = 0;
    int cavity_counter = 0;

    int processed = 0;
    size_t step_size = (birth_cells_.size() > 10) ? (birth_cells_.size() / 10) : 1;  // Show ~10 progress updates
    for (size_t i = 0; i < birth_cells_.size(); ++i) {
        if (i % step_size == 0) {
            std::cout << "    Processed " << i << " pairs for filtering..." << std::endl;
        }
        
        int birth_cell = birth_cells_[i];
        int death_cell = death_cells_[i];
        
        // Skip infinite persistence pairs
        if (death_cell == -1) {
            processed++;
            continue;
        }
        
        double birth_alpha = alphas[birth_cell];
        double death_alpha = alphas[death_cell];

        // Check if one alpha is positive and one is negative
        if ((birth_alpha!=death_alpha )&& (birth_alpha <= 0 && death_alpha > 0) && (birth_alpha >= core || death_alpha <= neighborhood)) {

            //also if birth cell is 0, skip
            if (birth_cell == 0) {
                processed++;
                std::cout << "Skipping 0 feature: " << birth_cell << " (alpha=" << birth_alpha << ")" << std::endl;
                std::cout << "Death cell: " << death_cell << " (alpha=" << death_alpha << ")" << std::endl;
                continue;
            }

            if (cell_dimensions[birth_cell] == 2) {
                cavity_counter++;
            }
            if (cell_dimensions[birth_cell] == 1) {
                handle_counter++;
            }
            if (cell_dimensions[birth_cell] == 0) {
                component_counter++;
            }
            // Get the generator directly from PHAT reduced matrix to avoid building Eigen matrix
            std::vector<int> generator;
            phat::column col_data;
            phat_matrix_.get_col(death_cell, col_data);
            generator.assign(col_data.begin(), col_data.end());
            std::sort(generator.begin(), generator.end());
            filtered_pairs_.push_back({birth_cell, death_cell, birth_alpha, death_alpha, generator});
        }

        processed++;
    }
    std::cout << "Betti numbers: " << component_counter << " / " << handle_counter << " / " << cavity_counter << std::endl;
    std::cout << "  Filtering completed. Found " << filtered_pairs_.size() << " filtered pairs" << std::endl;
}

void Persistence::convertFromPHAT() {
    try {
        std::cout << "  Starting PHAT to Eigen conversion..." << std::endl;
        int num_rows = boundary_matrix_.rows();
        int num_cols = boundary_matrix_.cols();
        std::cout << "  Matrix dimensions: " << num_rows << " x " << num_cols << std::endl;
        
        // Pre-allocate triplets vector with better size estimate
        std::vector<Eigen::Triplet<int>> triplets;
        // Estimate based on actual non-zero elements in PHAT matrix
        size_t estimated_nonzeros = 0;
        for (int col = 0; col < num_cols; ++col) {
            phat::column col_data;
            phat_matrix_.get_col(col, col_data);
            estimated_nonzeros += col_data.size();
        }
        triplets.reserve(estimated_nonzeros);
        
        // Convert to reduced matrix
        std::cout << "  Converting reduced matrix..." << std::endl;
        size_t step_size = (num_cols > 10) ? (num_cols / 10) : 1;  // Show ~10 progress updates
        for (int col = 0; col < num_cols; ++col) {
            if (col % step_size == 0) {
                std::cout << "    Processed " << col << " columns..." << std::endl;
            }
            phat::column col_data;
            phat_matrix_.get_col(col, col_data);
            for (int row : col_data) {
                if (row >= 0 && row < num_rows) {
                    triplets.emplace_back(row, col, 1);
                }
            }
        }

        reduced_matrix_.setFromTriplets(triplets.begin(), triplets.end());
        std::cout << "  Reduced matrix conversion completed. Nonzeros: " << reduced_matrix_.nonZeros() << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in convertFromPHAT: " << e.what() << std::endl;
        throw;
    }
}

const std::vector<PersistencePair>& Persistence::getFilteredPairs() const {
    return filtered_pairs_;
}

const std::vector<int>& Persistence::getBirthCells() const {
    return birth_cells_;
}

const std::vector<int>& Persistence::getDeathCells() const {
    return death_cells_;
}

const Eigen::SparseMatrix<int>& Persistence::getReducedMatrix() const {
    return reduced_matrix_;
}

const Eigen::SparseMatrix<int>& Persistence::getDualReducedMatrix() const {
    return dual_reduced_matrix_;
}

const phat::boundary_matrix<phat::bit_tree_pivot_column>& Persistence::getDualPHATMatrix() const {
    return dual_phat_matrix_;
}

std::vector<double> Persistence::outputPersistenceDiagram(const std::string& filename, double minX, double maxX, double minY, double maxY) const {
    // Extract base filename (before .wl)
    std::string base_filename = filename;
    size_t dot_pos = base_filename.find_last_of(".");
    if (dot_pos != std::string::npos && base_filename.substr(dot_pos) == ".wl") {
        base_filename = base_filename.substr(0, dot_pos);
    }
    
    const std::vector<double> alphas = boundary_matrix_.getAlphas();
    const std::vector<int>& cell_dimensions = boundary_matrix_.getDimensions();
    
    // First pass: collect all points by dimension and compute overall min/max
    std::vector<std::vector<std::pair<double, double>>> red_points_by_dim(3);   // birth <= 0 and death > 0
    std::vector<std::vector<std::pair<double, double>>> pink_points_by_dim(3);   // all others
    
    double overall_min_birth = 0.0, overall_max_birth = 0.0;
    double overall_min_death = 0.0, overall_max_death = 0.0;
    bool has_any_points = false;
    
    for (size_t i = 0; i < birth_cells_.size(); ++i) {
        int birth_cell = birth_cells_[i];
        int death_cell = death_cells_[i];
        
        // Check if birth_cell is valid
        if (birth_cell < 0 || birth_cell >= static_cast<int>(cell_dimensions.size())) {
            continue;
        }
        
        if (birth_cell >= static_cast<int>(alphas.size())) {
            continue;
        }
        
        // Get dimension of birth cell
        int dim = cell_dimensions[birth_cell];
        if (dim < 0 || dim > 2) {
            continue;  // Skip dimensions outside 0-2
        }
        
        double birth_alpha = alphas[birth_cell];
        double death_alpha;
        
        if (death_cell == -1) {
            // Skip infinite persistence pairs
            continue;
        } else if (death_cell >= 0 && death_cell < static_cast<int>(alphas.size())) {
            death_alpha = alphas[death_cell];
        } else {
            continue;
        }

        if (birth_alpha == death_alpha) {
            continue;
        }
        
        // Track overall min/max for plot range
        if (!has_any_points) {
            overall_min_birth = overall_max_birth = birth_alpha;
            overall_min_death = overall_max_death = death_alpha;
            has_any_points = true;
        } else {
            if (birth_alpha < overall_min_birth) overall_min_birth = birth_alpha;
            if (birth_alpha > overall_max_birth) overall_max_birth = birth_alpha;
            if (death_alpha < overall_min_death) overall_min_death = death_alpha;
            if (death_alpha > overall_max_death) overall_max_death = death_alpha;
        }
        
        // Categorize points by color and store by dimension
        if (birth_alpha <= 0 && death_alpha > 0) {
            red_points_by_dim[dim].push_back({birth_alpha, death_alpha});
        } else {
            pink_points_by_dim[dim].push_back({birth_alpha, death_alpha});
        }
    }
    // Calculate overall padding
    double overall_birth_padding = (overall_max_birth - overall_min_birth) * 0.1;
    double overall_death_padding = (overall_max_death - overall_min_death) * 0.1;
    if (overall_birth_padding == 0.0) overall_birth_padding = 0.1;
    if (overall_death_padding == 0.0) overall_death_padding = 0.1;


    double min_birth = overall_min_birth - overall_birth_padding;
    double max_birth = overall_max_birth + overall_birth_padding;
    double min_death = overall_min_death - overall_death_padding;
    double max_death = overall_max_death + overall_death_padding;

    if ((minX!=0) || (maxX!=0) || (minY!=0) || (maxY!=0)) {
        min_birth = std::min(min_birth, minX);
        max_birth = std::max(max_birth, maxX);
        min_death = std::min(min_death, minY);
        max_death = std::max(max_death, maxY);
    }

    double min_overall = std::min(min_birth, min_death);
    double max_overall = std::max(max_birth, max_death);
    // Second pass: output diagrams for each dimension using overall range
    for (int dim = 0; dim <= 2; ++dim) {
        std::string dim_filename = base_filename + "_dim" + std::to_string(dim) + ".wl";
        std::ofstream out_file(dim_filename);
        if (!out_file.is_open()) {
            std::cerr << "Error: Could not open file for writing: " << dim_filename << std::endl;
            continue;
        }
        
        const auto& red_points = red_points_by_dim[dim];
        const auto& pink_points = pink_points_by_dim[dim];
        bool has_points = !red_points.empty() || !pink_points.empty();
        
        // Write Mathematica code to create the persistence diagram
        out_file << "(* Persistence Diagram - Dimension " << dim << " *)\n";
        out_file << "(* Generated from persistence computation *)\n";
        out_file << "(* Contains all persistence pairs for dimension " << dim << " *)\n";
        out_file << "(* Red points: birth <= 0 and death > 0, Gray points: all others *)\n\n";
        
        // Write pink points
        out_file << "(* Gray points: all other points *)\n";
        out_file << "pinkPoints = {\n";
        for (size_t i = 0; i < pink_points.size(); ++i) {
            if (i > 0) out_file << ",\n";
            out_file << "  {" << pink_points[i].first << ", " << pink_points[i].second << "}";
        }
        out_file << "\n};\n\n";

        // Write red points
        out_file << "(* Red points: birth <= 0 and death > 0 *)\n";
        out_file << "redPoints = {\n";
        for (size_t i = 0; i < red_points.size(); ++i) {
            if (i > 0) out_file << ",\n";
            out_file << "  {" << red_points[i].first << ", " << red_points[i].second << "}";
        }
        out_file << "\n};\n\n";
        
        
        // Create the plot with lines for x=0 and y=0 using overall range
        out_file << "(* Create persistence diagram plot with x=0 and y=0 lines *)\n";
        if (has_points) {
            out_file << "persistenceDiagram = Show[\n";
            out_file << "  {\n";
            
            // Plot pink points
            if (!pink_points.empty()) {
                out_file << "    ListPlot[pinkPoints,\n";
                out_file << "      PlotStyle -> {PointSize[0.005], Gray}\n";
                out_file << "    ]";
            }

            // Plot red points
            if (!red_points.empty()) {
                if (!pink_points.empty()) {
                    out_file << ",\n";
                }
                out_file << "    ListPlot[redPoints,\n";
                out_file << "      PlotStyle -> {PointSize[0.01], Red}\n";
                out_file << "    ]";
            }
            
            out_file << "\n  },\n";
            out_file << "  PlotRange -> {{" << min_overall << ", " << max_overall << "}, {" << min_overall<< ", " << max_overall << "}},\n";
            out_file << "  AspectRatio -> 1,\n";
            out_file << "  Axes -> False,\n";
            out_file << " GridLines -> None,\n";
            out_file << "  Frame -> True,\n";
            out_file << "  FrameLabel -> {\"Birth\", \"Death\"},\n";
            out_file << "  PlotLabel -> \"Persistence Diagram - Dimension " << dim << "\",\n";
            out_file << "  Epilog -> {\n";

            out_file << "    {Black, Thickness[0.003], Arrow[{{" << min_overall << ", " << min_overall << "}, {"<< min_overall << ", " <<  max_death << "}}]},\n";
            out_file << "    {Black, Thickness[0.003], Arrow[{{" << min_overall << ", " << min_overall << "}, {"<< max_birth << ", " <<  min_overall << "}}]},\n";

            out_file << "    (* Line x = 0 *)\n";
            out_file << "    {Black, Thickness[0.003], Dashing[{0.01, 0.01}], Line[{{0, 0}, {0, " << max_death << "}}]},\n";
            out_file << "    (* Line y = 0 *)\n";
            out_file << "    {Black, Thickness[0.003], Dashing[{0.01, 0.01}], Line[{{0, 0}, {" << min_birth << ", 0}}]},\n";

            out_file << "    {GrayLevel[0.4], Thickness[0.002], Line[{{" << min_overall << ", " << min_overall << "}, {"<< max_overall << ", " <<  max_overall << "}}]}\n";

            out_file << "  }\n";
            out_file << "];\n\n";
        } else {
            // No points to plot
            out_file << "persistenceDiagram = Graphics[];\n\n";
        }
        
        out_file << "(* Display the plot *)\n";
        out_file << "persistenceDiagram\n";
        
        out_file.close();
        std::cout << "  Persistence diagram (dimension " << dim << ") written to: " << dim_filename << std::endl;
    }
    return {min_birth, max_birth, min_death, max_death};
} 