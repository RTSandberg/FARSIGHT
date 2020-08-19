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

#include "AMRStructure.hpp"

AMRStructure::AMRStructure() {}
AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0, 
                            int initial_height, 
                            double x_min, double x_max, double v_min, double v_max, 
                            double greens_epsilon, int num_steps, double dt, 
                            bool do_adaptively_refine)
                           : f0(f0), q(-1.0), qm(-1.0), 
                           initial_height(initial_height) , height(initial_height), max_height(initial_height),
                           x_min(x_min), x_max(x_max),
                           v_min(v_min), v_max(v_max), 
                           iter_num(0), num_steps(num_steps), dt(dt),
                           greens_epsilon(greens_epsilon), quad(trap),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                            do_adaptively_refine(do_adaptively_refine)
{
    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Lv = v_max - v_min;
    npanels_x = int(pow(2, initial_height));
    npanels_v = int(pow(2, initial_height));
    // create_prerefined_mesh();
    bool is_initial_step = true;
    
    generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
    
}
AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0, 
                            double q, double m, 
                            int initial_height, int max_height, 
                            double x_min, double x_max, double v_min, double v_max, 
                            double epsilon, Quadrature quad, int num_steps, double dt, 
                            bool do_adaptively_refine)
                           : f0(f0), q(q), qm(q/m), 
                           initial_height(initial_height), height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), v_min(v_min), v_max(v_max), 
                           iter_num(0), num_steps(num_steps), dt(dt),
                           greens_epsilon(greens_epsilon), quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine)
{
    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Lv = v_max - v_min;
    npanels_x = int(pow(2, initial_height));
    npanels_v = int(pow(2, initial_height));

    // create_prerefined_mesh();
    bool is_initial_step = true;
    generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
}
//end constructors

// getters
std::string AMRStructure::get_sim_dir() const { return sim_dir; }
// end getters

// setters

// end setters
