
#include <math.h> // M_PI, exp, cos, sin
#include <stdio.h> // printf
#include <string>
#include <string.h> // atof
#include <iostream> // cout, endl
#include <omp.h>
using std::cout;
using std::endl;
#include <stdexcept> // invalid_argument exception

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
    cout << "Hello world!" << endl;


    std::string sim_dir = "";

    double x_min = 0.0, x_max = 4 * M_PI, v_min = -6.0, v_max = 6.0;
    double Lx = x_max - x_min;
    double kx = 2.0 * M_PI / Lx;
    double amp = 0.5;
    double vth = 1.0;
    F0_LD f0{vth, kx, amp};
    
    int initial_height = 6; 
    double greens_epsilon = 0.2;
    int num_steps = 120;
    double dt = 0.5;
    bool do_adaptively_refine = false;

    AMRStructure amr{"/Users/ryansand/Documents/PlasmaProjects/Lagrangian_particle_method/heat_lamps_projects/amr/08_biquadratic_work/5_biquadratic_cpp/", f0, 
                initial_height, 
                x_min, x_max, v_min, v_max, 
                greens_epsilon, num_steps, dt, 
                do_adaptively_refine};

    amr.init_e();
    amr.write_to_file();

    for (int ii = 0; ii < num_steps; ++ii) {
        amr.step(false);
        amr.remesh();
        amr.write_to_file();
    }

    // cout << amr << endl;

}