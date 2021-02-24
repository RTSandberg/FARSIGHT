
#include <algorithm> // std::copy
#include <exception> // std::exception
#include <fstream>              // std::ifstream
#include <iterator> // std::ostream_iterator
#include <math.h> // M_PI, exp, cos, sin
#include <numeric> // std::inner_product
#include <stdio.h> // printf
#include <string>
#include <string.h> // atof
#include <iostream> // cout, endl
#include <omp.h>
using std::cout;
using std::endl;
#include <stdexcept> // invalid_argument exception
#include <thread> // std::thread::hardware_concurrency


extern "C" {
    #include <mpi.h>
}

#include "boost/property_tree/ptree.hpp"        //property_tree::ptree, property_tree::read_json
#include "boost/property_tree/json_parser.hpp"
namespace pt = boost::property_tree;

#include "AMRSimulation.hpp"
// #include "Panel.hpp"
// #include "AMRStructure.hpp"
// #include "initial_distributions.hpp"

// #define DEBUG // for debugging purposes
// #define DEBUG2


int main(int argc, char** argv) {

    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    std::string sim_dir, input_deck;
    if (argc > 1) {
        sim_dir = argv[1];
    } else {
        sim_dir = "";
    }
    if (argc > 3) {
        input_deck = argv[2];
        input_deck = input_deck + argv[3];
    } else {
        if (argc > 2) {
            input_deck = sim_dir + argv[2];
        } else {
            input_deck = sim_dir + "deck.json";
        }
    }

    AMRSimulation sim(sim_dir, input_deck);
    // initial gather
    // sim.step();
    sim.run();
    // evaluate e
    // write to diagnostic

    // sim.run();
    // step:
    // gather points
    // evaulate e, push 4 times
    // if remeshing:
    //     remesh
    //     initial gather and evaluate e
    // if not remeshing:
    //     evaluate e
    // if dumping to file:
    //     send field to points
    //     dump

    
    MPI_Finalize();
}
