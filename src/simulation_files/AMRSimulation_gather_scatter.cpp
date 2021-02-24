#include "AMRSimulation.hpp"

int AMRSimulation::gather() {
    
    xs = std::vector<double> ();
    vs = std::vector<double> ();
    q_ws = std::vector<double> ();
    es = std::vector<double> ();
    species_start = std::vector<size_t> ();
    species_end = std::vector<size_t> ();
    int Nx = 0;

    for (auto &species : species_list) {
        species_start.push_back(Nx);
        Nx += species->xs.size();
        species_end.push_back(Nx);
    }
    xs.reserve(Nx);
    vs.reserve(Nx);
    q_ws.reserve(Nx);
    es.reserve(Nx);
    
    for (auto &species : species_list) {
        xs.insert( xs.end(), species->xs.begin(), species->xs.end());
        vs.insert( vs.end(), species->vs.begin(), species->vs.end());
        q_ws.insert( q_ws.end(), species->q_ws.begin(), species->q_ws.end());
        es.insert( es.end(), species->es.begin(), species->es.end());
    }

    need_gather = false;
    return 0;
}

int AMRSimulation::scatter(bool send_e) {
    // distribute to each species
    size_t start_ind = 0;
    for (auto &species : species_list) {
        for (int ii = 0; ii < species->xs.size(); ++ii) {
            species->xs[ii] = xs[start_ind + ii];
            species->vs[ii] = vs[start_ind + ii];
            species->q_ws[ii] = q_ws[start_ind + ii];
            if (send_e) {
                species->es[ii] = es[start_ind + ii];
            }
        }
        start_ind += species->xs.size();
    }

    need_scatter = false;
    return 0;
}