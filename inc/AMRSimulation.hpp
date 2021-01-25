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

#include "boost/property_tree/ptree.hpp"        //property_tree::ptree, property_tree::read_json
#include "boost/property_tree/json_parser.hpp"
namespace pt = boost::property_tree;

// amr in-project dependencies
#include "initial_distributions.hpp"
#include "Panel.hpp"
#include "FieldStructure.hpp"
#include "AMRStructure.hpp"

struct AMRSimulation {
    std::string sim_dir;
    // domain parameters
    double Lx, Lv;
    double x_min, x_max;

    // time stepping parameters
    int iter_num;
    int num_steps;
    double dt;

    //field parameters
//     double greens_epsilon;
    Quadrature quad;
//     bool use_treecode;
    ElectricField* calculate_e;

    std::vector<AMRStructure> species_list;
    // ---- functions ---------
    AMRSimulation();
    int load_deck();
    int evaluate_field();
    int step();
    int run();
};
#endif /* AMRSIMULATION_HPP */
