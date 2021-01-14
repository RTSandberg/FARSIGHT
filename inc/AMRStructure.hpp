/**
 * @file AMRStructure.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-05-04
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#ifndef AMRSTRUCTURE_HPP
#define AMRSTRUCTURE_HPP

#include <algorithm>            // std::sort, std::copy, std::find
#include <assert.h>             // assert
#include <chrono>               // high_resolution _clock, duration_cast, microseconds
using namespace std::chrono;
#include <functional>           // std::logical_and
#include <fstream>
#include <iostream>             // std::cout, std::endl
#include <iterator>             // std::ostream_iterator for vector printing
#include <math.h>
#include <numeric>              // std::iota, std::accumulate
#include <omp.h>
#include <set>                  // std::set, set.find
#include <stdexcept>            // exceptions
#include <stdio.h>              // printf
#include <string> 
#include <vector>               // std::vector

// external libraries to include
#include <Eigen/Dense>
using namespace Eigen;
extern "C" {
    #include <mpi.h>
    #include <BaryTreeInterface.h>
}

// amr in-project dependencies
#include "initial_distributions.hpp"
#include "Panel.hpp"
#include "FieldStructure.hpp"


enum BoundaryConditions {periodic_bcs, open_bcs, last_bc};
enum Quadrature {trap, simpsons, last_quad};
enum ProfileTypes {sim_time, step_time, field_time, 
                remesh_time, tree_build_time, panel_test_time, interp_time, search_time, 
                eval_time, file_time, last_time};

struct AMRStructure {
    std::string sim_dir;
    // domain parameters
    double Lx, Lv;
    double x_min, x_max;
    double v_min, v_max;
    BoundaryConditions bcs;

    // species parameters
    double q, qm;

    // initial condition
    //std::function<double (double,double)> f0;
    distribution* f0;

    // mesh parameters
    int initial_height;
    int v_height;
    int height;
    int max_height;
    int npanels_x, npanels_v;
    double initial_dx, initial_dv;
    bool is_initial_mesh_set;
    bool need_further_refinement;
    int minimum_unrefined_index;
    bool do_adaptively_refine;
    std::vector<double> amr_epsilons;

    double Q0, M0;

    std::vector <Panel> panels;
    std::vector <int> leaf_inds;
    // std::vector <Particle> particles;
    std::vector<double> xs, vs, fs, q_ws, es;

    std::vector <Panel> old_panels;
    std::vector<double> old_xs, old_vs, old_fs;
    // class InterpolateDistribution : public Distribution {
    //     double operator() (double x, double v);
    // }
    // InterpolateDistribution interp_f;
    // interpolation parameters
    bool allow_boundary_extrapolation = false;
    // bool do_unshear = false;
    bool do_unshear = false;
    // bool sqrt_f = false;
    bool sqrt_f = false;

    // time stepping parameters
    int iter_num;
    int num_steps;
    double dt;

    //field parameters
//     double greens_epsilon;
    Quadrature quad;
//     bool use_treecode;
    ElectricField* calculate_e;

    //profile parameters
    bool do_profile;
    std::vector<duration<double>> time_operations;
    std::vector<int> num_operations;

    // private functions
    int create_prerefined_mesh();
    int create_prerefined_mesh_v_refinement();
    void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine);
    void refine_panels_refine_v(std::function<double (double,double)> f, bool do_adaptive_refine);
    void test_panel(int panel_ind, bool verbose);

    int write_particles_to_file();
    int write_panels_to_file();
    int write_particles_to_file(bool pre_remesh);
    int write_panels_to_file(bool pre_remesh);

    public:
        AMRStructure();
        AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)>& f0, 
                int initial_height, 
                double x_min, double x_max, double v_min, double v_max, 
                ElectricField* calculate_e, int num_steps, double dt);
        AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)>& f0, 
                int initial_height, int max_height,
                double x_min, double x_max, double v_min, double v_max, 
                BoundaryConditions bcs,
                ElectricField* calculate_e, int num_steps, double dt,
                bool do_adaptively_refine, std::vector<double>& amr_epsilons);
        AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)>& f0, 
                double q, double m, 
                int initial_height, int v_height, int max_height, 
                double x_min, double x_max, double v_min, double v_max, 
                BoundaryConditions bcs,
                ElectricField* calculate_e, Quadrature quad, int num_steps, double dt, 
                bool do_adaptively_refine, std::vector<double>& amr_epsilons);

        // end constructors
        // destructor
        ~AMRStructure();
        // getters
        std::vector<double> get_e();
        std::string get_sim_dir() const;
        // end getters

        // amr
        void generate_mesh(std::function<double (double,double)> f,
                        bool do_adaptive_refine, bool is_initial_step);
        void set_leaves_weights();
        void recursively_set_leaves_weights(int panel_ind);

        // remesh
        void copy_to_old();
        void reset_mesh();
        void remesh();

        // interpolation functions
        bool use_limiter;
        double limit_val;
        void shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& vs);
        int find_leaf_containing_xv_recursively(double &x, const double &v, bool& beyond_boundary, int panel_ind, bool verbose);
        int find_leaf_containing_point_from_neighbor(double& tx, double& tv, bool& beyond_boundary, int leaf_ind, std::set<int>& history, bool verbose);
        // int find_leaf_containing();
        void interpolate_to_initial_xvs(std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& vs, int nx, int nv,bool verbose);
        double interpolate_from_mesh(double xs, double vs, bool verbose);
        void interpolate_from_mesh(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        void interpolate_from_mesh_slow(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        double interpolate_from_panel(double x, double v, int panel_ind,bool verbose);
        void interpolate_from_panel_to_points(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& vs,
                                                std::vector<int>& point_inds, int panel_ind, bool use_limiter, double limit_val);

        // field functions
        // void calculate_e();
        // void vector_calc_e_wrapper();
        void init_e();

        // step functions
        void step(bool get_4th_e);

        // io
        friend std::ostream& operator<<(std::ostream& os, const AMRStructure& amr);
        void print_amr();
        int write_to_file();
        int write_to_file(bool pre_remesh);
        void print_panel_points();

        // profiling
        void add_time(ProfileTypes prof_type, duration<double> op_time);
        void print_times();
};



// int calculate_E_mq(double* es, const double* targets, int nt, 
//                   const double* sources, const double* q_ws, int ns, 
//                   double L, double epsilon);

#endif /* AMRSTRUCTURE_HPP */

