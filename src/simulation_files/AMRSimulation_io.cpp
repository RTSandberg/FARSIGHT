#include "AMRSimulation.hpp"

void AMRSimulation::print_sim_setup() {

    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << "deck found in: " << deck_address << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << v_min << " <= v <= " << v_max << endl;
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
    field_object->print_field_obj();
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