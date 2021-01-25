#include "AMRSimulation.hpp"

int AMRSimulation::run() {

    auto start = high_resolution_clock::now();
    amr.init_e();
    auto stop = high_resolution_clock::now();
    amr.add_time(field_time, duration_cast<duration<double>>(stop - start) );
    amr.write_to_file();

    for (int ii = 1; ii < num_steps+1; ++ii) {
        bool get_4th_e = false;
        start = high_resolution_clock::now();
        amr.step(get_4th_e);
        stop = high_resolution_clock::now();
        amr.add_time(step_time, duration_cast<duration<double>>(stop - start) );


        if ((ii) % n_steps_remesh == 0) {

#ifdef DEBUG
bool pre_remesh = true;
amr.write_to_file(pre_remesh);
#endif

            start = high_resolution_clock::now();
            amr.remesh();
            stop = high_resolution_clock::now();
            amr.add_time(remesh_time, duration_cast<duration<double>>(stop - start) );
        }


        if ((ii) % n_steps_diag == 0) {

            auto file_start = high_resolution_clock::now();
            amr.write_to_file();
            auto file_stop = high_resolution_clock::now();
            amr.add_time(file_time,  duration_cast<duration<double>>(file_stop - file_start) );
        }
    //     // cout << amr << endl;
    }

    // cout << amr << endl;
    auto sim_stop = high_resolution_clock::now();
    amr.add_time(sim_time,  duration_cast<duration<double>>(sim_stop - sim_start) );
    // cout << "Sim time " << sim_duration.count() << " seconds" << endl;

    cout << "End with " << amr.panels.size() << " panels, " << amr.xs.size() << " points." << endl; 

    return 0;
    }
