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
    double Lx, Lv;
    double x_min, x_max, v_min, v_max;

    // time stepping parameters
    int iter_num;
    int num_steps;
    int n_steps_remesh;
    int n_steps_diag;
    double dt;

    //field parameters
    // double greens_epsilon;
    BoundaryConditions bcs;
    Quadrature quad;
    bool use_treecode;
    ElectricField* field_object;

    std::vector<distribution*> ic_list;
    std::vector<AMRStructure*> species_list;
    // ---- functions ---------
    // constructor
    AMRSimulation();
    // destructor
    ~AMRSimulation();
    // ------------------
    AMRSimulation(std::string sim_dir, std::string deck_address);
    int load_deck(std::string &deck_address, pt::ptree &deck);
    int get_box_t_params(pt::ptree &deck);
    distribution* make_f0_return_ptr(pt::ptree &species_deck_portion);
    ElectricField* make_field_return_ptr(pt::ptree &deck);
    AMRStructure* make_species_return_ptr(pt::ptree &species_deck_portion, distribution* f0);
    int evaluate_field_uniform_grid();
    int evaluate_field();
    int step();
    int run();
    void print_sim_setup();
};
#endif /* AMRSIMULATION_HPP */
