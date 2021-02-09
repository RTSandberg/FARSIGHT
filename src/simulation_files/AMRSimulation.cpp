#include "AMRSimulation.hpp"

AMRSimulation::AMRSimulation() {}

AMRSimulation::AMRSimulation(std::string sim_dir, std::string deck_address) {

    sim_dir = sim_dir;
    deck_address = deck_address;
    // cout << "deck location: " << deck_address << endl;

    // create ptree
    pt::ptree deck;
    load_deck(deck_address, deck);

    get_box_t_params(deck);

    // create e solver
    // load species

    // 
    // create species list
    

    // print AMR description
    print_sim_setup();
}