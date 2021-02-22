#include "AMRSimulation.hpp"

int AMRSimulation::step() {
    iter_num += 1;
    std::cout << "step " << iter_num << std::endl;



    if (need_gather) {
        // gather from species (need at first and after every remesh)
        gather();
        need_gather = false;
    }

    // rk4 step
    rk4_step(false);

    need_scatter = true;

    // if remesh: remesh_and_calculate_e, scatter if needed, set need_gather to true
    // if not remesh: evaluate e (remeshing ends by calculating e on uniform grid)
    // if dump : scatter if needed, write to file

    return 0;
}

int AMRSimulation::rk4_step(bool get_4th_e) {



// initialize rk4 vectors
    std::vector<double> xtemp = xs;
    std::vector<double> v1 = vs, v2, v3, v4;
    std::vector<double> a1 = es, a2(xs.size()), a3(xs.size()), a4(xs.size());
    int N = xs.size();


//   math_vector v1, v2, v3, v4;
//   math_vector f1(total_num_points), f2(total_num_points), f3(total_num_points), f4(total_num_points);
//   math_vector tempx;
  
    // k1 = (v1,a1) = F(un) = F(xn,pn) = (vn, q/m E(xn) )
  // v1 = vs = vn
    // calculate_E_mq(a1, xs, xs, q_ws, L, epsilon);
    // v1 = vs;
    for (size_t sp_i; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            a1[xi] *= species_qms[sp_i];
        }
    }
    
    // k2 = (v2,a2) = F(un + h/2 k1) = F(xn + delt/2 v1, vn + delt/2 a1)
    //              = ( vn + delt a1 / 2, q/m E(xn + delt v1 /2) )
    for (int ii = 0; ii < N; ++ii) {
        v2.push_back(vs[ii] + 0.5 * dt * a1[ii]);
        xtemp[ii] += 0.5 * dt * v1[ii];
    }

    // auto start = high_resolution_clock::now();
    
    std::vector<double> xtemp_cpy (xtemp);
    std::vector<double> xtemp_cpy2 (xtemp);
    (*calculate_e)(a2.data(), xtemp_cpy.data(), a2.size(),
                    xtemp_cpy2.data(), q_ws.data(), xtemp.size());
    // auto stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );
    for (size_t sp_i; sp_i < N_sp; species_list) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            a2[xi] *= species_qms[sp_i];
        }
    }

    // k3 = (v3,a3) = F(un + h/2 k2) = F(xn + delt/2 v2, vn + delt/2 a2)
    //              = ( vn + delt a2 / 2, q/m E(xn + delt v2 /2) )
    for (int ii = 0; ii < N; ++ii) {
        v3.push_back(vs[ii] + 0.5 * dt * a2[ii]);
        xtemp[ii] = xs[ii] + 0.5 * dt * v2[ii];
    }
    // start = high_resolution_clock::now();
    xtemp_cpy = xtemp;
    xtemp_cpy2 = xtemp;
    (*calculate_e)(a3.data(), xtemp_cpy.data(), a3.size(),
                    xtemp_cpy2.data(), q_ws.data(), xtemp.size());
    // stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );
    // for (int ii = 0; ii < N; ++ii) {
    //     a3[ii] *= qm;
    // }
    for (size_t sp_i; sp_i < N_sp; species_list) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            a3[xi] *= species_qms[sp_i];
        }
    }
    // k4 = (v4,a4) = F(un + h k3) = F(xn + delt v3, pn + delt f3)
    //              = ( (pn + delt f3) / m, q E(xn + delt v3) )
    for (int ii = 0; ii < N; ++ii) {
        v4.push_back(vs[ii] + dt * a3[ii]);
        xtemp[ii] = xs[ii] + dt * v3[ii];
    }
    // start = high_resolution_clock::now();

    xtemp_cpy = xtemp;
    xtemp_cpy2 = xtemp;
    (*calculate_e)(a4.data(), xtemp_cpy.data(), a4.size(), 
                    xtemp_cpy2.data(), q_ws.data(), xtemp.size());
    // stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );

    // for (int ii = 0; ii < N; ++ii) {
    //     a4[ii] *= qm;
    // }
    for (size_t sp_i; sp_i < N_sp; species_list) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            a4[xi] *= species_qms[sp_i];
        }
    }
    // un+1 = (xn+1,pn+1) = un + h/6 (k1 + 2k2 + 2k3 + k4)
    //                    = (xn + delt/6 (v1 + 2 v2 + 2 v3 + v4),
    //                    = pn + delt/6 (a1 + 2 a2 + 2 a3 + a4) )
    // store xn+1, pn+1 in xn,pn
    for (int ii = 0; ii < N; ++ii) {
        xs[ii] += dt / 6.0 * (v1[ii] + 2 * v2[ii] + 2 * v3[ii] + v4[ii]);
        vs[ii] += dt / 6.0 * (a1[ii] + 2 * a2[ii] + 2 * a3[ii] + a4[ii]);
    }
    if (get_4th_e) {
        // start = high_resolution_clock::now();
        xtemp_cpy = xs;
        xtemp_cpy2 = xs;
        (*calculate_e)(es.data(), xtemp_cpy.data(), xs.size(),
                        xtemp_cpy2.data(), q_ws.data(), xs.size());
        // stop = high_resolution_clock::now();
        // add_time(field_time, duration_cast<duration<double>>(stop - start) );
    }
    

    return 0;
}