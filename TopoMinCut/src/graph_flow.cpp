#include "graph_flow.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <queue>
#include <limits>

std::pair<std::vector<int>, std::vector<int>>
compute_min_cut(
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
) {
    using namespace boost;
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string,
            property < vertex_index_t, long,
                property < vertex_color_t, boost::default_color_type,
                    property < vertex_distance_t, long,
                        property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        property < edge_capacity_t, float,
            property < edge_residual_capacity_t, float,
                property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

    Graph g;
    auto capacity = get(edge_capacity, g);
    auto residual_capacity = get(edge_residual_capacity, g);
    auto rev = get(edge_reverse, g);

    Traits::vertex_descriptor s = add_vertex(g);
    std::vector<Traits::vertex_descriptor> birth_nodes(real_birth_cells.size());
    for (size_t i = 0; i < real_birth_cells.size(); ++i)
        birth_nodes[i] = add_vertex(g);
    std::vector<Traits::vertex_descriptor> death_nodes(real_death_cells.size());
    for (size_t j = 0; j < real_death_cells.size(); ++j)
        death_nodes[j] = add_vertex(g);
    Traits::vertex_descriptor t = add_vertex(g);

    float total_sum = 0.0f;
    for (int idx : real_birth_cells) total_sum += std::abs(alphas[idx]);
    for (int idx : real_death_cells) total_sum += std::abs(alphas[idx]);

    // Section 1: s -> birth nodes
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        float cap = total_sum - std::abs(alphas[real_birth_cells[i]]);
        auto e = add_edge(s, birth_nodes[i], g).first;
        auto rev_e = add_edge(birth_nodes[i], s, g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;
    }

    // Section 2: birth nodes -> death nodes (only if conflict)
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        for (size_t j = 0; j < real_death_cells.size(); ++j) {
            if (conflicts[i][j] == 1) {
                auto e = add_edge(birth_nodes[i], death_nodes[j], g).first;
                auto rev_e = add_edge(death_nodes[j], birth_nodes[i], g).first;
                capacity[e] = std::numeric_limits<float>::max();
                capacity[rev_e] = 0;
                residual_capacity[e] = std::numeric_limits<float>::max();
                residual_capacity[rev_e] = 0;
                rev[e] = rev_e;
                rev[rev_e] = e;
            }
        }
    }

    // Section 3: death nodes -> t
    for (size_t j = 0; j < real_death_cells.size(); ++j) {
        float cap = total_sum - std::abs(alphas[real_death_cells[j]]);
        auto e = add_edge(death_nodes[j], t, g).first;
        auto rev_e = add_edge(t, death_nodes[j], g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;
    }

    // Compute max flow
    boykov_kolmogorov_max_flow(g, s, t);

    // BFS to find min-cut
    std::vector<bool> visited(num_vertices(g), false);
    std::queue<Traits::vertex_descriptor> q;
    q.push(s);
    visited[s] = true;

    while (!q.empty()) {
        auto u = q.front(); q.pop();
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            if (residual_capacity[e] > 0 && !visited[target(e, g)]) {
                visited[target(e, g)] = true;
                q.push(target(e, g));
            }
        }
    }

    std::vector<int> cut_birth_nodes, cut_death_nodes;
    for (auto u : boost::make_iterator_range(vertices(g))) {
        if (!visited[u]) continue;
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            auto v = target(e, g);
            if (visited[u] && !visited[v] && capacity[e] > 0) {
                for (size_t i = 0; i < birth_nodes.size(); ++i) {
                    if (u == s && v == birth_nodes[i]) {
                        cut_birth_nodes.push_back(i);
                    }
                }
                for (size_t j = 0; j < death_nodes.size(); ++j) {
                    if (u == death_nodes[j] && v == t) {
                        cut_death_nodes.push_back(j);
                    }
                }
            }
        }
    }
    return {cut_birth_nodes, cut_death_nodes};
} 



std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_with_dimensions(
    const std::vector<int>& dimensions,
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
) {
    using namespace boost;
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string,
            property < vertex_index_t, long,
                property < vertex_color_t, boost::default_color_type,
                    property < vertex_distance_t, long,
                        property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        property < edge_capacity_t, float,
            property < edge_residual_capacity_t, float,
                property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

    Graph g;
    auto capacity = get(edge_capacity, g);
    auto residual_capacity = get(edge_residual_capacity, g);
    auto rev = get(edge_reverse, g);

    Traits::vertex_descriptor s = add_vertex(g);
    std::vector<Traits::vertex_descriptor> birth_nodes(real_birth_cells.size());
    for (size_t i = 0; i < real_birth_cells.size(); ++i)
        birth_nodes[i] = add_vertex(g);
    std::vector<Traits::vertex_descriptor> death_nodes(real_death_cells.size());
    for (size_t j = 0; j < real_death_cells.size(); ++j)
        death_nodes[j] = add_vertex(g);
    Traits::vertex_descriptor t = add_vertex(g);

    float total_sum = 0.0f;

    float total_sum_used_for_handle_cost = 0.0f;

    for (int idx : real_birth_cells){

        total_sum_used_for_handle_cost += std::abs(alphas[idx]);
        
    }
    for (int idx : real_death_cells) total_sum_used_for_handle_cost += std::abs(alphas[idx]);

    total_sum_used_for_handle_cost += 1.0f;


    for (int idx : real_birth_cells){
        if (dimensions[idx] == 1) {
            total_sum += total_sum_used_for_handle_cost;
            //total_sum += std::abs(alphas[idx]);
        } else {
            total_sum += std::abs(alphas[idx]);
        }
    }
    for (int idx : real_death_cells){
        if (dimensions[idx] == 1) {
            //total_sum += std::abs(alphas[idx]);
            total_sum += std::abs(alphas[idx]);
        } else {
            total_sum += std::abs(alphas[idx]);
        }
    }

    // Section 1: s -> birth nodes
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        float cap = 0;
        if (dimensions[real_birth_cells[i]] == 1) {
            cap = total_sum - total_sum_used_for_handle_cost;
            //cap = total_sum - std::abs(alphas[real_birth_cells[i]]);
        } else {
            cap = total_sum - std::abs(alphas[real_birth_cells[i]]);
        }
        auto e = add_edge(s, birth_nodes[i], g).first;
        auto rev_e = add_edge(birth_nodes[i], s, g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;
    }

    // Section 2: birth nodes -> death nodes (only if conflict)
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        for (size_t j = 0; j < real_death_cells.size(); ++j) {
            if (conflicts[i][j] == 1) {
                auto e = add_edge(birth_nodes[i], death_nodes[j], g).first;
                auto rev_e = add_edge(death_nodes[j], birth_nodes[i], g).first;
                capacity[e] = std::numeric_limits<float>::max();
                capacity[rev_e] = 0;
                residual_capacity[e] = std::numeric_limits<float>::max();
                residual_capacity[rev_e] = 0;
                rev[e] = rev_e;
                rev[rev_e] = e;
            }
        }
    }

    // Section 3: death nodes -> t
    for (size_t j = 0; j < real_death_cells.size(); ++j) {
        float cap = 0;
        if (dimensions[real_death_cells[j]] == 1 ) {
            //cap = total_sum - total_sum_used_for_handle_cost;
            cap = total_sum - std::abs(alphas[real_death_cells[j]]);
        } else {
            cap = total_sum - std::abs(alphas[real_death_cells[j]]);
        }
        auto e = add_edge(death_nodes[j], t, g).first;
        auto rev_e = add_edge(t, death_nodes[j], g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;
    }

    // Compute max flow
    boykov_kolmogorov_max_flow(g, s, t);

    // BFS to find min-cut
    std::vector<bool> visited(num_vertices(g), false);
    std::queue<Traits::vertex_descriptor> q;
    q.push(s);
    visited[s] = true;

    while (!q.empty()) {
        auto u = q.front(); q.pop();
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            if (residual_capacity[e] > 0 && !visited[target(e, g)]) {
                visited[target(e, g)] = true;
                q.push(target(e, g));
            }
        }
    }

    std::vector<int> cut_birth_nodes, cut_death_nodes;
    for (auto u : boost::make_iterator_range(vertices(g))) {
        if (!visited[u]) continue;
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            auto v = target(e, g);
            if (visited[u] && !visited[v] && capacity[e] > 0) {
                for (size_t i = 0; i < birth_nodes.size(); ++i) {
                    if (u == s && v == birth_nodes[i]) {
                        cut_birth_nodes.push_back(i);
                    }
                }
                for (size_t j = 0; j < death_nodes.size(); ++j) {
                    if (u == death_nodes[j] && v == t) {
                        cut_death_nodes.push_back(j);
                    }
                }
            }
        }
    }
    return {cut_birth_nodes, cut_death_nodes};
} 


std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_simple_cost_function(
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
) {
    using namespace boost;
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string,
            property < vertex_index_t, long,
                property < vertex_color_t, boost::default_color_type,
                    property < vertex_distance_t, long,
                        property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        property < edge_capacity_t, float,
            property < edge_residual_capacity_t, float,
                property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

    Graph g;
    auto capacity = get(edge_capacity, g);
    auto residual_capacity = get(edge_residual_capacity, g);
    auto rev = get(edge_reverse, g);

    Traits::vertex_descriptor s = add_vertex(g);
    std::vector<Traits::vertex_descriptor> birth_nodes(real_birth_cells.size());
    for (size_t i = 0; i < real_birth_cells.size(); ++i)
        birth_nodes[i] = add_vertex(g);
    std::vector<Traits::vertex_descriptor> death_nodes(real_death_cells.size());
    for (size_t j = 0; j < real_death_cells.size(); ++j)
        death_nodes[j] = add_vertex(g);
    Traits::vertex_descriptor t = add_vertex(g);

    float total_sum = 0.0f;

    
    for (size_t i = 0; i < real_birth_cells.size(); ++i){
        total_sum += real_generating_sets_fore_indices[i].size();
    } 
    for (size_t j = 0; j < real_death_cells.size(); ++j){
        total_sum += real_generating_sets_back_indices[j].size();
    }

    // Section 1: s -> birth nodes
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        float cap = total_sum - real_generating_sets_fore_indices[i].size();
        auto e = add_edge(s, birth_nodes[i], g).first;
        auto rev_e = add_edge(birth_nodes[i], s, g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;

    }

    // Section 2: birth nodes -> death nodes (only if conflict)
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        for (size_t j = 0; j < real_death_cells.size(); ++j) {
            if (conflicts[i][j] == 1) {
                auto e = add_edge(birth_nodes[i], death_nodes[j], g).first;
                auto rev_e = add_edge(death_nodes[j], birth_nodes[i], g).first;
                capacity[e] = std::numeric_limits<float>::max();
                capacity[rev_e] = 0;
                residual_capacity[e] = std::numeric_limits<float>::max();
                residual_capacity[rev_e] = 0;
                rev[e] = rev_e;
                rev[rev_e] = e;
            }
        }
    }

    // Section 3: death nodes -> t
    for (size_t j = 0; j < real_death_cells.size(); ++j) {
        float cap = total_sum - real_generating_sets_back_indices[j].size();
        auto e = add_edge(death_nodes[j], t, g).first;
        auto rev_e = add_edge(t, death_nodes[j], g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;

    }

    // Compute max flow
    boykov_kolmogorov_max_flow(g, s, t);

    // BFS to find min-cut
    std::vector<bool> visited(num_vertices(g), false);
    std::queue<Traits::vertex_descriptor> q;
    q.push(s);
    visited[s] = true;

    while (!q.empty()) {
        auto u = q.front(); q.pop();
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            if (residual_capacity[e] > 0 && !visited[target(e, g)]) {
                visited[target(e, g)] = true;
                q.push(target(e, g));
            }
        }
    }

    std::vector<int> cut_birth_nodes, cut_death_nodes;
    for (auto u : boost::make_iterator_range(vertices(g))) {
        if (!visited[u]) continue;
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            auto v = target(e, g);
            if (visited[u] && !visited[v] && capacity[e] > 0) {
                for (size_t i = 0; i < birth_nodes.size(); ++i) {
                    if (u == s && v == birth_nodes[i]) {
                        cut_birth_nodes.push_back(i);
                    }
                }
                for (size_t j = 0; j < death_nodes.size(); ++j) {
                    if (u == death_nodes[j] && v == t) {
                        cut_death_nodes.push_back(j);
                    }
                }
            }
        }
    }
    return {cut_birth_nodes, cut_death_nodes};
} 


std::pair<std::vector<int>, std::vector<int>>
compute_min_cut_simple_cost_function_with_dimensions(
    const std::vector<int>& dimensions,
    const std::vector<int>& real_birth_cells,
    const std::vector<int>& real_death_cells,
    const std::vector<std::vector<int>>& real_generating_sets_fore_indices,
    const std::vector<std::vector<int>>& real_generating_sets_back_indices,
    const std::vector<std::vector<int>>& conflicts,
    const std::vector<double>& alphas
) {
    using namespace boost;
    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string,
            property < vertex_index_t, long,
                property < vertex_color_t, boost::default_color_type,
                    property < vertex_distance_t, long,
                        property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        property < edge_capacity_t, float,
            property < edge_residual_capacity_t, float,
                property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

    Graph g;
    auto capacity = get(edge_capacity, g);
    auto residual_capacity = get(edge_residual_capacity, g);
    auto rev = get(edge_reverse, g);

    Traits::vertex_descriptor s = add_vertex(g);
    std::vector<Traits::vertex_descriptor> birth_nodes(real_birth_cells.size());
    for (size_t i = 0; i < real_birth_cells.size(); ++i)
        birth_nodes[i] = add_vertex(g);
    std::vector<Traits::vertex_descriptor> death_nodes(real_death_cells.size());
    for (size_t j = 0; j < real_death_cells.size(); ++j)
        death_nodes[j] = add_vertex(g);
    Traits::vertex_descriptor t = add_vertex(g);

    float total_sum = 0.0f;
    float total_sum_used_for_handle_cost = 0.0f;

        
    for (size_t i = 0; i < real_birth_cells.size(); ++i){
        if (dimensions[real_birth_cells[i]] == 1) {
            total_sum_used_for_handle_cost += real_generating_sets_fore_indices[i].size();
        } else {
            total_sum_used_for_handle_cost += real_generating_sets_fore_indices[i].size();
        }
    } 
    for (size_t j = 0; j < real_death_cells.size(); ++j){
        total_sum_used_for_handle_cost += real_generating_sets_back_indices[j].size();
    }

    total_sum_used_for_handle_cost += 1.0f;


    for (size_t i = 0; i < real_birth_cells.size(); ++i){
        if (dimensions[real_birth_cells[i]] == 1) {
            total_sum += real_generating_sets_fore_indices[i].size();
        } else {
            total_sum += real_generating_sets_fore_indices[i].size();
        }
    } 
    for (size_t j = 0; j < real_death_cells.size(); ++j){
        total_sum += real_generating_sets_back_indices[j].size();
    }

    // Section 1: s -> birth nodes
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {

        float cap = 0;
        if (dimensions[real_birth_cells[i]] == 1) {
            cap = total_sum - real_generating_sets_fore_indices[i].size();
        } else {
            cap = total_sum - real_generating_sets_fore_indices[i].size();
        }
        auto e = add_edge(s, birth_nodes[i], g).first;
        auto rev_e = add_edge(birth_nodes[i], s, g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;

    }

    // Section 2: birth nodes -> death nodes (only if conflict)
    for (size_t i = 0; i < real_birth_cells.size(); ++i) {
        for (size_t j = 0; j < real_death_cells.size(); ++j) {
            if (conflicts[i][j] == 1) {
                auto e = add_edge(birth_nodes[i], death_nodes[j], g).first;
                auto rev_e = add_edge(death_nodes[j], birth_nodes[i], g).first;
                capacity[e] = std::numeric_limits<float>::max();
                capacity[rev_e] = 0;
                residual_capacity[e] = std::numeric_limits<float>::max();
                residual_capacity[rev_e] = 0;
                rev[e] = rev_e;
                rev[rev_e] = e;
            }
        }
    }

    // Section 3: death nodes -> t
    for (size_t j = 0; j < real_death_cells.size(); ++j) {
        float cap = total_sum - real_generating_sets_back_indices[j].size();
        auto e = add_edge(death_nodes[j], t, g).first;
        auto rev_e = add_edge(t, death_nodes[j], g).first;
        capacity[e] = cap;
        capacity[rev_e] = 0;
        residual_capacity[e] = cap;
        residual_capacity[rev_e] = 0;
        rev[e] = rev_e;
        rev[rev_e] = e;

    }

    // Compute max flow
    boykov_kolmogorov_max_flow(g, s, t);

    // BFS to find min-cut
    std::vector<bool> visited(num_vertices(g), false);
    std::queue<Traits::vertex_descriptor> q;
    q.push(s);
    visited[s] = true;

    while (!q.empty()) {
        auto u = q.front(); q.pop();
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            if (residual_capacity[e] > 0 && !visited[target(e, g)]) {
                visited[target(e, g)] = true;
                q.push(target(e, g));
            }
        }
    }

    std::vector<int> cut_birth_nodes, cut_death_nodes;
    for (auto u : boost::make_iterator_range(vertices(g))) {
        if (!visited[u]) continue;
        for (auto e : boost::make_iterator_range(out_edges(u, g))) {
            auto v = target(e, g);
            if (visited[u] && !visited[v] && capacity[e] > 0) {
                for (size_t i = 0; i < birth_nodes.size(); ++i) {
                    if (u == s && v == birth_nodes[i]) {
                        cut_birth_nodes.push_back(i);
                    }
                }
                for (size_t j = 0; j < death_nodes.size(); ++j) {
                    if (u == death_nodes[j] && v == t) {
                        cut_death_nodes.push_back(j);
                    }
                }
            }
        }
    }
    return {cut_birth_nodes, cut_death_nodes};
} 