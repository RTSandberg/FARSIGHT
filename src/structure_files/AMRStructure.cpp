/**
 * @file AMRStructure.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2020-05-05
 * 
 * @copyright Copyright (c) 2020
 * 
 */

/*
New features
---
allow for p_levels, int, >= 0 
*/

#include "AMRStructure.hpp"

AMRStructure::AMRStructure() {}
AMRStructure::AMRStructure(std::string sim_dir,std::string species_name,
                            distribution* f0, //std::function<double (double,double)> f0, 
                            int initial_height, 
                            double x_min, double x_max, double p_min, double p_max)
                           : f0(f0), q(-1.0), qm(-1.0), 
                           initial_height(initial_height) , p_height(0),
                           height(initial_height), max_height(initial_height),
                           x_min(x_min), x_max(x_max),
                           p_min(p_min), p_max(p_max), bcs(periodic_bcs), quad(trap),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(false),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    this->species_name = species_name;
    Lx = x_max - x_min;
    Lp = p_max - p_min;
    npanels_x = int(pow(2, initial_height));
    npanels_p = int(pow(2, initial_height + p_height));
    // create_prerefined_mesh();
    bool is_initial_step = true;
    
    generate_mesh([&](double x, double p) { return (*f0)(x,p); }, do_adaptively_refine, is_initial_step);
}

        // AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)>& f0, 
        //         int initial_height, 
        //         double x_min, double x_max, double p_min, double p_max, 
        //         ElectricField* calculate_e, int num_steps, double dt,
        //         bool do_adaptively_refine, std::vector<double>& amr_epsilons);
AMRStructure::AMRStructure(std::string sim_dir, std::string species_name,
                            distribution* f0, //std::function<double (double,double)> f0,
                            int initial_height, int max_height, 
                            double x_min, double x_max, double p_min, double p_max, 
                            BoundaryConditions bcs,
                            bool do_adaptively_refine, std::vector<double>& amr_epsilons)
                           : f0(f0), q(-1.0), qm(-1.0), 
                           initial_height(initial_height), p_height(0),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), p_min(p_min), p_max(p_max), bcs(bcs),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    this->species_name = species_name;
    Lx = x_max - x_min;
    Lp = p_max - p_min;
    npanels_x = int(pow(2, initial_height));
    npanels_p = int(pow(2, initial_height + p_height));
    initial_dx = Lx / npanels_x;
    initial_dp = Lp / npanels_p;
    this->amr_epsilons = amr_epsilons;

    // create_prerefined_mesh();
    bool is_initial_step = true;
    generate_mesh([&](double x, double p) { return (*f0)(x,p); }, do_adaptively_refine, is_initial_step);
    // calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
}

AMRStructure::AMRStructure(std::string sim_dir, std::string species_name,
                            distribution* f0, //std::function<double (double,double)> f0, 
                            double q, double m, 
                            int initial_height, int p_height, int max_height, 
                            double x_min, double x_max, double p_min, double p_max, 
                            BoundaryConditions bcs, Quadrature quad,
                            bool do_adaptively_refine, std::vector<double>& amr_epsilons)
                           : f0(f0), q(q), qm(q/m), 
                           initial_height(initial_height), p_height(p_height),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), p_min(p_min), p_max(p_max), bcs(bcs),
                           quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    this->species_name = species_name;
    Lx = x_max - x_min;
    Lp = p_max - p_min;
    npanels_x = int(pow(2, initial_height));
    npanels_p = int(pow(2, initial_height + p_height));
    initial_dx = Lx / npanels_x;
    initial_dp = Lp / npanels_p;
    this->amr_epsilons = amr_epsilons;

    // create_prerefined_mesh();
    bool is_initial_step = true;
    generate_mesh([&](double x, double p) { return (*f0)(x,p); }, do_adaptively_refine, is_initial_step);
    // calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
}
//end constructors

//destructor
AMRStructure::~AMRStructure() = default;

// getters
// std::vector<double> AMRStructure::get_e() { return es; };
std::string AMRStructure::get_sim_dir() const { return sim_dir; }
// end getters

// setters

// end setters

void AMRStructure::add_time(ProfileTypes prof_type, duration<double> op_time) {
    num_operations[prof_type] ++;
    time_operations[prof_type] += op_time;
}


void AMRStructure::get_reduced_xs_ws() {
    std::vector<size_t> inds (xs.size());

    std::iota(inds.begin(), inds.end(), 0);
    std::sort(inds.begin(), inds.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });

    // unique pat
    int s_ind = inds[0];
    double val = xs[s_ind];
    reduced_xs = std::vector<double> (1,val);
    inv_inds_reduced_xs = std::vector<size_t> (xs.size());
    int unique_ind = 0;
    for (int ii = 1; ii < xs.size(); ++ii) {
        s_ind = inds[ii];
        double xsi = xs[s_ind];
        // std::cout << "ii " << ii << ", x[sort_ind[ii]] " << xsi << ", cf " << val << std::endl;
        // if (xsi > val) { std::cout << "xsi > val";}
        // else {std::cout << "xsi <= val" << std::endl;}
        if (xsi > val) {
            val = xsi;
            unique_ind++;
            reduced_xs.push_back(xsi);
        } 
        // std::cout << "inv_inds[" << s_ind << "] set to " << unique_ind << std::endl;
        inv_inds_reduced_xs[s_ind] = unique_ind;
    }

    // std::vector<double> sort_ws (species_reduced_xs.size());
    reduced_ws = std::vector<double> (reduced_xs.size());
    for (int ii = 0; ii < inv_inds_reduced_xs.size(); ++ii) {
        reduced_ws[inv_inds_reduced_xs[ii]] += q_ws[ii];
    }
}

void AMRStructure::get_reduced_es(double* reduced_es) {
    es = std::vector<double>(xs.size());
    for (int ii = 0; ii < inv_inds_reduced_xs.size(); ++ii) {
        // es.push_back(sort_es[inv_inds[ii]]);
        es[ii] = reduced_es[inv_inds_reduced_xs[ii]];
    }

}


