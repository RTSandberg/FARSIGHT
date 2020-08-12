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

#include <algorithm> // sort, copy, std::find
#include <assert.h> // assert
#include <chrono> // high_resolution _clock, duration_cast, microseconds
using namespace std::chrono;
#include <functional>
#include <fstream>
#include <iostream> // std::cout, std::endl
#include <iterator> // ostream_iterator for vector printing
#include <math.h>
#include <numeric> // iota
#include <omp.h>
#include <set> // std::set, set.find
#include <stdexcept> // exceptions
#include <stdio.h> // printf
#include <string> 
#include <vector> // std::vector

#include <Eigen/Dense>
using namespace Eigen;

#include "Panel.hpp"

enum Quadrature {simpsons, trap};

struct AMRStructure {
    std::string sim_dir;
    // domain parameters
    double Lx, Lv;
    double x_min, x_max;
    double v_min, v_max;

    // species parameters
    double q, qm;

    // initial condition
    std::function<double (double,double)> f0;

    // mesh parameters
    int initial_height;
    int height;
    int max_height;
    int npanels_x, npanels_v;
    double initial_dx, initial_dv;
    bool is_initial_mesh_set;
    bool need_further_refinement;
    int minimum_unrefined_index;
    bool do_adaptively_refine;

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

    // time stepping parameters
    int iter_num;
    int num_steps;
    double dt;

    //field parameters
    double greens_epsilon;
    Quadrature quad;

    // private functions
    int create_prerefined_mesh();
    void refine_panels(std::function<double (double,double)> f, bool do_adaptive_refine);
    void test_panel(int panel_ind);

    int write_particles_to_file();
    int write_panels_to_file();

    public:
        AMRStructure();
        AMRStructure(std::string sim_dir, std::function<double (double,double)> f0, 
                int initial_height, 
                double x_min, double x_max, double v_min, double v_max, 
                double greens_epsilon, int num_steps, double dt, 
                bool do_adaptively_refine);
        AMRStructure(std::string sim_dir, std::function<double (double,double)> f0, 
                double q, double m, 
                int initial_height, int max_height, 
                double x_min, double x_max, double v_min, double v_max, 
                double greens_epsilon, Quadrature quad, int num_steps, double dt, 
                bool do_adaptively_refine);

        // end constructors
        // getters
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
        void shift_xs(std::vector<double>& shifted_xs, const std::vector<double>& xs, const std::vector<double>& vs);
        int find_leaf_containing_xv_recursively(double &x, const double &v, int panel_ind, bool verbose);
        int find_leaf_containing_point_from_neighbor(double& tx, double& tv, int leaf_ind, std::set<int>& history, bool verbose);
        // int find_leaf_containing();
        void interpolate_to_initial_xvs(std::vector<double>& fs, std::vector<double>& xs, std::vector<double>& vs, int nx, int nv,bool verbose);
        double interpolate_from_mesh(double xs, double vs, bool verbose);
        void interpolate_from_mesh(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        void interpolate_from_mesh_slow(std::vector<double> &values, std::vector<double>& x, std::vector<double>& v, bool verbose);
        double interpolate_from_panel(double x, double v, int panel_ind,bool verbose);
        void interpolate_from_panel_to_points(std::vector<double>& values, std::vector<double>& xs, std::vector<double>& vs,
                                                std::vector<int>& point_inds, int panel_ind);

        // field functions
        void calculate_e();
        void vector_calc_e_wrapper();
        void init_e();

        // step functions
        void step(bool get_4th_e);

        // io
        friend std::ostream& operator<<(std::ostream& os, const AMRStructure& amr);
        void print_amr();
        int write_to_file();
        void print_panel_points();

};



int calculate_E_mq(std::vector<double>& es, const std::vector<double>& targets, const std::vector<double>& sources, const std::vector<double>& q_weights, double L, double epsilon);

#endif /* AMRSTRUCTURE_HPP */

