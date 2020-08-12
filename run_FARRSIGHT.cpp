
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

class F0_LD {
    double vth, k, amp;
    public:
        F0_LD() : vth(1.0), k(.5), amp(.1) {}
        F0_LD(double vth, double k, double amp) : vth(vth), k(k), amp(amp) {}

        double operator() (double x, double v) {
            return 1.0 / sqrt(2.0 * M_PI) / vth * exp(-v*v / 2 / vth /vth) * ( 1 + amp * cos(k * x ));
        }
        double get_vth() { return vth; }
        double get_k() { return k; }
        double get_amp() { return amp; }
};

int main(int argc, char** argv) {

    if (argc < 14) {
        std::cout << "Not enough input arguments: had " << argc - 1 << ", but need 13" << endl;
        std::cout << "Usage: 1:sim_dir 2:xmin 3:xmax 4:vmin 5:vmax 6:normal_k 7:amp 8:vth";
        std::cout << " 9:initial_height 10:greens_epsilon";
        std::cout << " 11:num_steps 12:n_steps_remesh 13:dt" << std::endl;
        return 1;
    }
    std::string sim_dir = argv[1];
    // sim_dir.append("/");

    double x_min = atof(argv[2]), x_max = atof(argv[3]), v_min = atof(argv[4]), v_max = atof(argv[5]);
    double Lx = x_max - x_min;
    double kx = 2.0 * M_PI / Lx * atof(argv[6]);
    double amp = atof(argv[7]);//0.5;
    double vth = atof(argv[8]);//1.0;
    F0_LD f0{vth, kx, amp};
    
    int initial_height = atoi(argv[9]);//6; 
    double greens_epsilon = atof(argv[10]);//0.2;
    int num_steps = atoi(argv[11]);//120;
    int n_steps_remesh = atoi(argv[12]);
    double dt = atof(argv[13]);//0.5;
    bool do_adaptively_refine = false;


    cout << "============================" << endl;
    cout << "Running a FRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << v_min << " <= v <= " << v_max << endl;
    cout << "k=" << kx << ", amp = " << amp << ", vth = " << vth << endl;
    cout << "height " << initial_height << endl;
    cout << "green's epsilon = " << greens_epsilon << endl;

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

        cout << "Step time " << duration.count() << " microseconds." << endl << endl;

        start = high_resolution_clock::now();
        amr.remesh();
        stop = high_resolution_clock::now();
        duration = duration_cast<microseconds>(stop - start);
        cout << "remesh time " << duration.count() << " microseconds." << endl << endl;

    // cout << "Old data copy time " << duration.count() << " microseconds." << endl << endl;
        amr.write_to_file();
    //     // cout << amr << endl;
    }

    // cout << amr << endl;

}