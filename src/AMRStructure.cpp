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
allow for v_levels, int, >= 0 
*/

#include "AMRStructure.hpp"

AMRStructure::AMRStructure() {}
AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0, 
                            int initial_height, 
                            double x_min, double x_max, double v_min, double v_max, 
                            ElectricField* calculate_e, int num_steps, double dt)
                           : f0(f0), q(-1.0), qm(-1.0), 
                           initial_height(initial_height) , v_height(0),
                           height(initial_height), max_height(initial_height),
                           x_min(x_min), x_max(x_max),
                           v_min(v_min), v_max(v_max), bcs(periodic_bcs),
                           iter_num(0), num_steps(num_steps), dt(dt),
                           calculate_e(calculate_e), quad(trap),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(false),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Lv = v_max - v_min;
    npanels_x = int(pow(2, initial_height));
    npanels_v = int(pow(2, initial_height + v_height));
    // create_prerefined_mesh();
    bool is_initial_step = true;
    
    generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
    f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
    cout << "extrapolating value is " << f_beyond_boundary << endl;
}

        // AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)>& f0, 
        //         int initial_height, 
        //         double x_min, double x_max, double v_min, double v_max, 
        //         ElectricField* calculate_e, int num_steps, double dt,
        //         bool do_adaptively_refine, std::vector<double>& amr_epsilons);
AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0,
                            int initial_height, int max_height, 
                            double x_min, double x_max, double v_min, double v_max, 
                            BoundaryConditions bcs,
                            ElectricField* calculate_e, int num_steps, double dt, 
                            bool do_adaptively_refine, std::vector<double>& amr_epsilons)
                           : f0(f0), q(-1.0), qm(-1.0), 
                           initial_height(initial_height), v_height(0),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), v_min(v_min), v_max(v_max), bcs(bcs),
                           iter_num(0), num_steps(num_steps), dt(dt),
                           calculate_e(calculate_e),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Lv = v_max - v_min;
    npanels_x = int(pow(2, initial_height));
    npanels_v = int(pow(2, initial_height + v_height));
    initial_dx = Lx / npanels_x;
    initial_dv = Lv / npanels_v;
    this->amr_epsilons = amr_epsilons;

    // create_prerefined_mesh();
    bool is_initial_step = true;
    generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
    // calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
    f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
    cout << "extrapolating value is " << f_beyond_boundary << endl;
}

AMRStructure::AMRStructure(std::string sim_dir, distribution* f0, //std::function<double (double,double)> f0, 
                            double q, double m, 
                            int initial_height, int v_height, int max_height, 
                            double x_min, double x_max, double v_min, double v_max, 
                            BoundaryConditions bcs,
                            ElectricField* calculate_e, Quadrature quad, int num_steps, double dt, 
                            bool do_adaptively_refine, std::vector<double>& amr_epsilons)
                           : f0(f0), q(q), qm(q/m), 
                           initial_height(initial_height), v_height(v_height),
                           height(initial_height), max_height(max_height), 
                           x_min(x_min), x_max(x_max), v_min(v_min), v_max(v_max), bcs(bcs),
                           iter_num(0), num_steps(num_steps), dt(dt),
                           calculate_e(calculate_e), quad(quad),
                           is_initial_mesh_set(false), minimum_unrefined_index(0), need_further_refinement(false),
                           do_adaptively_refine(do_adaptively_refine),
                           use_limiter(false), limit_val(0.0)
{
    time_operations = std::vector<duration<double>>(last_time);
    num_operations = std::vector<int> (last_time);

    this->sim_dir = sim_dir;
    Lx = x_max - x_min;
    Lv = v_max - v_min;
    npanels_x = int(pow(2, initial_height));
    npanels_v = int(pow(2, initial_height + v_height));
    initial_dx = Lx / npanels_x;
    initial_dv = Lv / npanels_v;
    this->amr_epsilons = amr_epsilons;

    // create_prerefined_mesh();
    bool is_initial_step = true;
    generate_mesh([&](double x, double v) { return (*f0)(x,v); }, do_adaptively_refine, is_initial_step);
    // calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
    f_beyond_boundary = *std::min_element(fs.begin(), fs.end() );
    cout << "extrapolating value is " << f_beyond_boundary << endl;
}
//end constructors

//destructor
AMRStructure::~AMRStructure() = default;

// getters
std::vector<double> AMRStructure::get_e() { return es; };
std::string AMRStructure::get_sim_dir() const { return sim_dir; }
// end getters

// setters

// end setters

void AMRStructure::add_time(ProfileTypes prof_type, duration<double> op_time) {
    num_operations[prof_type] ++;
    time_operations[prof_type] += op_time;
}
