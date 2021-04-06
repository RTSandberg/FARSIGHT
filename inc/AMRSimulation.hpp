/**
 * @file AMRSimulation.hpp
 * @author Ryan Sandberg (ryan.t.sandberg@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-01-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef AMRSIMULATION_HPP
#define AMRSIMULATION_HPP

#include <algorithm>            // std::sort, std::copy, std::find
#include <assert.h>             // assert
#include <chrono>               // high_resolution _clock, duration_cast, microseconds
using namespace std::chrono;
#include <functional>           // std::logical_and
#include <fstream>
#include <iostream>             // std::cout, std::endl
using std::cout;
using std::endl;
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
// #include <Eigen/Dense>
// using namespace Eigen;
extern "C" {
    #include <mpi.h>
    #include <BaryTreeInterface.h>
}

#include "boost/property_tree/ptree.hpp"        //property_tree::ptree, property_tree::read_json
#include "boost/property_tree/json_parser.hpp"
namespace pt = boost::property_tree;

// amr in-project dependencies
// #include "initial_distributions.hpp"
// #include "Panel.hpp"
#include "FieldStructure.hpp"
#include "AMRStructure.hpp"


// enum BoundaryConditions {periodic_bcs, open_bcs, last_bc};
// // enum Quadrature {trap, simpsons, last_quad};
// enum ProfileTypes {sim_time, step_time, field_time, 
//                 remesh_time, tree_build_time, panel_test_time, interp_time, search_time, 
//                 eval_time, file_time, amr_test_time, amr_refine_time, last_time};

struct AMRSimulation {
    std::string sim_dir;
    std::string deck_address;
    // domain parameters
    double Lx, Lp;
    double x_min, x_max, p_min, p_max;

    // time stepping parameters
    bool relativistic;
    int iter_num;
    double t;
    int num_steps;
    int n_steps_remesh;
    int n_steps_diag;
    double dt;
    bool need_scatter;
    bool need_gather;

    //field parameters
    // double greens_epsilon;
    BoundaryConditions bcs;
    Quadrature quad;
    bool use_treecode;
    ElectricField* calculate_e;
    bool use_external_field;
    ExternalElectricField* calculate_e_external;

    std::vector<distribution*> ic_list;
    std::vector<AMRStructure*> species_list;
    int N_sp;
    std::vector<size_t> species_start, species_end;

    std::vector<double> xs, ps, q_ws, es;
    std::vector<double> species_qs, species_ms, species_qms;
    // ---- functions ---------
    // constructor
    AMRSimulation();
    AMRSimulation(std::string sim_dir, std::string deck_address);
    // destructor
    ~AMRSimulation();
    // ------------------
    int load_deck(std::string &deck_address, pt::ptree &deck);
    int get_box_t_params(pt::ptree &deck);
    distribution* make_f0_return_ptr(pt::ptree &species_deck_portion);
    ElectricField* make_field_return_ptr(pt::ptree &deck);
    void make_external_field(pt::ptree &deck);
    AMRStructure* make_species_return_ptr(pt::ptree &species_deck_portion, distribution* f0);
    void get_qms();
    int gather();
    int scatter(bool send_e);

    int evaluate_field_uniform_grid(double t);
    int evaluate_field(std::vector<double>& es_local, std::vector<double>& xs_local, std::vector<double>& q_ws_local, double t);


    void relativistic_momentum_push(double* xtemp, double* vtemp, double* ptemp);
    void nonrelativistic_momentum_push(double* xtemp, double* vtemp, double* ptemp);
    int step();
    int rk4_step(bool get_4th_e);

    int remesh();
    int run();
    int write_to_file();
    int write_to_file(bool pre_remesh);
    void print_sim_setup();
};
#endif /* AMRSIMULATION_HPP */
