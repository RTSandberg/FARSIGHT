#include "AMRSimulation.hpp"

int AMRSimulation::write_to_file() { 
    bool pre_remesh = false;
    return write_to_file(pre_remesh);
}
int AMRSimulation::write_to_file(bool pre_remesh) {

    if (need_scatter) {
        bool send_e = true;
        scatter(send_e);
    }
    for (auto& species : species_list) {
        species->write_to_file(pre_remesh, iter_num);
    }
    return 0;
}

void AMRSimulation::print_sim_setup() {

    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << "deck found in: " << deck_address << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << p_min << " <= p <= " << p_max << endl;
    switch (bcs) {
        case (open_bcs) : cout << "Using open boundary conditions" << endl;
            break;
        default : // periodic
            cout << "Using periodic boundary conditions" << endl;
            break;
    }
    // cout << "height " << initial_height << ", v height " << v_height << endl;
    switch (quad) {
        case (simpsons) : cout << "Using Simpson's rule" << endl;
            break;
        default : cout << "Using trap rule" << endl;
            break;
    }
    cout << "Taking " << num_steps << " steps with dt = " << dt << endl;
    cout << "Remesh every " << n_steps_remesh << " step(s), diagnostic dump every " << n_steps_diag << " step(s)" << endl;

// print field data
    calculate_e->print_field_obj();
    if (use_external_field) {
        calculate_e_external->print_field_obj();
    }
    // cout << "green's epsilon = " << greens_epsilon << endl;
    // cout << "use treecode flag " << use_treecode << endl;
    // if (use_treecode > 0) { 
    //     if (0 <= beta && beta <= 1.0) {
    //         cout << "Using treecode with beta " << beta << endl;
    //     } else {
    //         cout << "Using treecode with mac " << mac << " and degree " << degree << endl;
    //     }
    // } else {
    //     cout << "using direct sum" << endl;
    // }
    // print species data
    // }
    for (int ii = 0; ii < species_list.size(); ++ii) {
        cout << *(species_list[ii]) << endl;
    }
    cout << "============================" << endl;

}