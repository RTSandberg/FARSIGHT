
#include <algorithm> // std::copy
#include <iterator> // std::ostream_iterator
#include <math.h> // M_PI, exp, cos, sin
#include <numeric> // std::inner_product
#include <stdio.h> // printf
#include <string>
#include <string.h> // atof
#include <iostream> // cout, endl
// #include <mpi.h>
#include <omp.h>
using std::cout;
using std::endl;
#include <stdexcept> // invalid_argument exception
#include <thread> // std::thread::hardware_concurrency

extern "C" {
    #include <mpi.h>
}

#include "Panel.hpp"
#include "AMRStructure.hpp"
#include "initial_distributions.hpp"




int main(int argc, char** argv) {

    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (argc < 23) {
        std::cout << "Not enough input arguments: had " << argc - 1 << ", but need 13" << std::endl;
        std::cout << "Usage: 1:sim_dir 2:xmin 3:xmax 4:vmin 5:vmax" << std::endl;
        std::cout << " 6:sim_type 7:normal_k 8:amp 9:vth 10:vstr" << std::endl;
        std::cout << " 11:initial_height 12:greens_epsilon" << std::endl;
        std::cout << " 13:use_treecode 14:treecode_beta" << std::endl;
        std::cout << " 15: mac 16: degree 17: max_source 18: max target" << std::endl;
        std::cout << " 19:num_steps 20:n_steps_remesh 21: n_steps_diag 22:dt" << std::endl;
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
    int sim_type = atoi(argv[6]);

    distribution* f0;

    switch (sim_type)
    {
        case 1: // weak Landau Damping
            f0 = new F0_LD(vth, kx, amp);
            break;
        case 2: // strong Landau Damping
            f0 = new F0_LD(vth, kx, amp);
            break;
        case 3: // 'strong' two-stream
            f0 = new F0_strong_two_stream(vth, kx, amp);
            break;
        case 4: // 'colder' two-stream
            f0 = new F0_colder_two_stream(vth, vstr, kx, amp);
            break;
        default:
            f0 = new F0_LD(vth, kx, amp);
            break;
    }

    // F0_colder_two_stream f0{vth, vstr, kx, amp};
    
    int initial_height = atoi(argv[11]);//6; 
    double greens_epsilon = atof(argv[12]);//0.2;
    int use_treecode = atoi(argv[13]);
    double beta = atof(argv[14]);
    double mac = atof(argv[15]);
    int degree = atoi(argv[16]);
    int max_source = atoi(argv[17]);
    int max_target = atoi(argv[18]);

    // int nxsqrt = pow(2, initial_height + 1) + 1;
    // int nx = nxsqrt * nxsqrt;

    ElectricField* calculate_e;
    if (use_treecode > 0) {
        if (0 <= beta && beta <= 1.0) {
            calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, beta);
        } else {
            int verbosity = 0;
            calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, 
                                            mac, degree, 
                                            max_source, max_target, 
                                            verbosity);
        }
    } else {
        calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
    }

    int num_steps = atoi(argv[19]);//120;
    int n_steps_remesh = atoi(argv[20]);
    int n_steps_diag = atoi(argv[21]);
    double dt = atof(argv[22]);//0.5;
    bool do_adaptively_refine = false;


    cout << "============================" << endl;
    cout << "Running a FARRSIGHT simulation" << endl;
    cout << "sim dir: " << sim_dir << endl;
    cout << x_min << " <= x <= " << x_max << endl;
    cout << v_min << " <= v <= " << v_max << endl;
    cout << "k=" << kx << ", amp = " << amp << ", vth = " << vth << ", vstr = " << vstr <<  endl;
    cout << "height " << initial_height << endl;
    cout << "green's epsilon = " << greens_epsilon << endl;
    cout << "use treecode flag " << use_treecode << endl;
    if (use_treecode > 0) { 
        if (0 <= beta && beta <= 1.0) {
            cout << "Using treecode with beta " << beta << endl;
        } else {
            cout << "Using treecode with mac " << mac << " and degree " << degree << endl;
        }
    } else {
        cout << "using direct sum" << endl;
    }
    cout << "============================" << endl;

    auto sim_start = high_resolution_clock::now();

    AMRStructure amr{sim_dir, f0, 
                initial_height, 
                x_min, x_max, v_min, v_max, 
                calculate_e, num_steps, dt, 
                do_adaptively_refine};

// ------ problem with tc!  
/*
    cout << "Problem with treecode.  Trying to debug" << endl;

    E_MQ_DirectSum ds {Lx, greens_epsilon};
    AMRStructure amr_ds{sim_dir, f0, 
                initial_height, 
                x_min, x_max, v_min, v_max, 
                &ds, num_steps, dt, 
                do_adaptively_refine};

    amr.init_e();
    amr_ds.init_e();

    bool get_4th_e = false;
    amr.step(get_4th_e);
    amr_ds.step(get_4th_e);
    amr.remesh();
    amr_ds.remesh();

    std::vector<double> es_ds = amr_ds.get_e();
    std::vector<double> es_tc = amr.get_e();

    int nt = es_ds.size();
    double dx_t = Lx / nt;
    std::vector<double> error(nt);
    double l2norm_ds=0.0, l2norm_tc=0.0, l2norm_diff=0.0;
    for (int ii = 0; ii < nt; ++ii) {
        double dsii = es_ds[ii];
        double tcii = es_tc[ii];
        l2norm_ds += dsii * dsii * dx_t;
        l2norm_tc += tcii * tcii * dx_t;
        double errii = es_ds[ii] - es_tc[ii];
        error[ii] = errii;
        l2norm_diff += errii * errii * dx_t;
    }
    // l2norm_ds = std::inner_product(es_ds.begin(), es_ds.end(), es_ds.begin(), 0);
    // l2norm_tc = std::inner_product(es_tc.begin(), es_tc.end(), es_tc.begin(), 0);
    // l2norm_diff = std::inner_product(error.begin(), error.end(), error.begin(), 0);
    l2norm_ds = sqrt(l2norm_ds);
    l2norm_tc = sqrt(l2norm_tc);
    l2norm_diff = sqrt(l2norm_diff);

    cout << "l2 norm of es, direct sum = " << l2norm_ds << endl;
    cout << "l2 norm of es, treecode = " << l2norm_tc << endl;
    cout << "l2 norm of error = " << l2norm_diff << endl;
*/
// ----- end treecode debug section


    auto start = high_resolution_clock::now();
    amr.init_e();
    auto stop = high_resolution_clock::now();
    amr.add_time(field_time, duration_cast<duration<double>>(stop - start) );
    amr.write_to_file();

    for (int ii = 0; ii < num_steps; ++ii) {
        bool get_4th_e = false;
        start = high_resolution_clock::now();
        amr.step(get_4th_e);
        stop = high_resolution_clock::now();
        amr.add_time(step_time, duration_cast<duration<double>>(stop - start) );


        start = high_resolution_clock::now();
        amr.remesh();
        stop = high_resolution_clock::now();
        amr.add_time(remesh_time, duration_cast<duration<double>>(stop - start) );

        if ((ii+1) % n_steps_diag == 0) {
            amr.write_to_file();
        }
    //     // cout << amr << endl;
    }

    // cout << amr << endl;
    auto sim_stop = high_resolution_clock::now();
    amr.add_time(sim_time,  duration_cast<duration<double>>(sim_stop - sim_start) );
    // cout << "Sim time " << sim_duration.count() << " seconds" << endl;

    amr.print_times();
    

    delete f0;
    delete calculate_e;
    
    MPI_Finalize();
}
