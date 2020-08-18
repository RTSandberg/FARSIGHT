
#include <math.h> // M_PI, exp, cos, sin
#include <stdio.h> // printf
#include <string>
#include <string.h> // atof
#include <iostream> // cout, endl
#include <omp.h>
using std::cout;
using std::endl;
#include <stdexcept> // invalid_argument exception
#include <thread> // std::thread::hardware_concurrency

#include "Panel.hpp"
#include "AMRStructure.hpp"
#include "initial_distributions.hpp"




int main(int argc, char** argv) {

    if (argc < 14) {
        std::cout << "Not enough input arguments: had " << argc - 1 << ", but need 13" << std::endl;
        std::cout << "Usage: 1:sim_dir 2:xmin 3:xmax 4:vmin 5:vmax " << std::endl;
        std::cout << "6:sim_type 7:normal_k 8:amp 9:vth 10:vstr" << std::endl;
        std::cout << " 11:initial_height 12:greens_epsilon" << std::endl;
        std::cout << " 13:num_steps 14:n_steps_remesh 15: n_steps_diag 16:dt" << std::endl;
        return 1;
    }
    std::string sim_dir = argv[1];
    // sim_dir.append("/");

    double x_min = atof(argv[2]), x_max = atof(argv[3]), v_min = atof(argv[4]), v_max = atof(argv[5]);
    double Lx = x_max - x_min;
    double kx = 2.0 * M_PI / Lx * atof(argv[7]);
    double amp = atof(argv[8]);//0.5;
    double vth = atof(argv[9]);//1.0;
    double vstr = atof(argv[10]);
    // switch (atof(argv[6]))
    // {
    // case 1:
        
    //     break;
    
    // default:
    //     break;
    // }
    F0_colder_two_stream f0{vth, vstr, kx, amp};
    
    int initial_height = atoi(argv[11]);//6; 
    double greens_epsilon = atof(argv[12]);//0.2;
    int num_steps = atoi(argv[13]);//120;
    int n_steps_remesh = atoi(argv[14]);
    int n_steps_diag = atoi(argv[15]);
    double dt = atof(argv[16]);//0.5;
    bool do_adaptively_refine = false;


    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << v_min << " <= v <= " << v_max << endl;
    cout << "k=" << kx << ", amp = " << amp << ", vth = " << vth << ", vstr = " << vstr <<  endl;
    cout << "height " << initial_height << endl;
    cout << "green's epsilon = " << greens_epsilon << endl;

    auto sim_start = high_resolution_clock::now();

    AMRStructure amr{sim_dir, f0, 
                initial_height, 
                x_min, x_max, v_min, v_max, 
                greens_epsilon, num_steps, dt, 
                do_adaptively_refine};

    amr.init_e();
    amr.write_to_file();

    // cout << amr << endl;


    // bool get_4th_e = true;
    // amr.step(get_4th_e);
    // cout << amr << endl;
    // cout << "Now trying to remesh" << endl;
    // amr.remesh();
    // amr.write_to_file();

    for (int ii = 0; ii < num_steps; ++ii) {
        bool get_4th_e = false;
        auto start = high_resolution_clock::now();
        amr.step(get_4th_e);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        // cout << "Step time " << duration.count() << " microseconds." << endl << endl;

        start = high_resolution_clock::now();
        amr.remesh();
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        // cout << "remesh time " << duration.count() << " microseconds." << endl << endl;

    // cout << "Old data copy time " << duration.count() << " microseconds." << endl << endl;
        amr.write_to_file();
    //     // cout << amr << endl;
    }

    // cout << amr << endl;
    auto sim_stop = high_resolution_clock::now();
    auto sim_duration = duration_cast<microseconds>(sim_stop - sim_start);
    cout << "Sim time " << sim_duration.count() << endl;

}
