#include "AMRSimulation.hpp"

int AMRSimulation::remesh() {
    if (need_scatter) {
        //scatter to species
        bool send_e = false;
        scatter(send_e);
    }

    // remesh
    for (auto &species : species_list) {
        species->remesh();
    }

    evaluate_field_uniform_grid(t);

    need_gather = true;
    return 0;
}