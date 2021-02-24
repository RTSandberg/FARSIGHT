#include "AMRSimulation.hpp"


int AMRSimulation::evaluate_field_uniform_grid() {
    // get reduced xs, reduced_ws from each species
    // combine reduced x, reduced_ws
    std::vector<double> reduced_xs, reduced_ws;
    for (auto &species : species_list) {
        species->get_reduced_xs_ws();
        reduced_xs.insert( reduced_xs.end(), species->reduced_xs.begin(), species->reduced_xs.end());
        reduced_ws.insert( reduced_ws.end(), species->reduced_ws.begin(), species->reduced_ws.end());
    }

    // duplicate reduced x
    std::vector<double> reduced_xs_cpy (reduced_xs);
    // evaluate reduced e
    std::vector<double> reduced_es(reduced_ws.size());
    (*calculate_e)(reduced_es.data(), reduced_xs.data(), reduced_xs.size(),
                 reduced_xs_cpy.data(), reduced_ws.data(), reduced_xs.size() );
    // distribute to each species
    size_t start_ind = 0;
    for (auto &species : species_list) {
        species->get_reduced_es(reduced_es.data() + start_ind);
        start_ind += species->reduced_xs.size();
    }
    need_gather = true;
    return 0;
}

int AMRSimulation::evaluate_field() {


    return 0;
}