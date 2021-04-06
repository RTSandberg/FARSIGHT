#include "AMRSimulation.hpp"

// #define DEBUG

int AMRSimulation::step() {

    if (need_gather) {
        // gather from species (need at first and after every remesh)
        gather();
    }

    // rk4 step
    rk4_step(false);

    iter_num += 1;
    t += dt;
    std::cout << "step " << iter_num << std::endl;

    // if remesh: remesh_and_calculate_e, scatter if needed, set need_gather to true
    if (iter_num % n_steps_remesh == 0) {
        bool save_Lagrangian = true;
        if (save_Lagrangian) {
            write_to_file(true);
        }
        remesh();
    }
    // if not remesh: evaluate e (remeshing ends by calculating e on uniform grid)
    else {
        evaluate_field(es, xs, q_ws, t);
        // std::vector<double> xtemp_cpy (xs);
        // std::vector<double> xtemp_cpy2 (xs);
        // (*calculate_e)(es.data(), xtemp_cpy.data(), es.size(),
        //                 xtemp_cpy2.data(), q_ws.data(), xtemp_cpy.size());
    }

    // if dump : write to file
    if (iter_num % n_steps_diag == 0) {
        write_to_file();
    }

    return 0;
}

void AMRSimulation::relativistic_momentum_push(double* xtemp, double* vtemp, double* ptemp) {
    for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            double norm_p = ptemp[xi]/ species_ms[sp_i];
            double gamma_1_i = sqrt(1 + norm_p * norm_p);
            double vi = norm_p / gamma_1_i;
            vtemp[xi] = vi;
            xtemp[xi] += 0.5 * dt * vi;
        }
    }
}
void AMRSimulation::nonrelativistic_momentum_push(double* xtemp, double* vtemp, double* ptemp) {
    for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            double vi = ptemp[xi]/ species_ms[sp_i];
            vtemp[xi] = vi;
            xtemp[xi] += 0.5 * dt * vi;
        }
    }
}

int AMRSimulation::rk4_step(bool get_4th_e) {

// initialize rk4 vectors
    std::vector<double> xk = xs;
    std::vector<double> v1(xs.size()), v2(xs.size()), v3(xs.size()), v4(xs.size() );
    std::vector<double> pk = ps;
    std::vector<double> ef1 = es, ef2(xs.size()), ef3(xs.size()), ef4(xs.size());
    int N = xs.size();


//   math_vector p1, p2, p3, p4;
//   math_vector f1(total_num_points), f2(total_num_points), f3(total_num_points), f4(total_num_points);
//   math_vector tempx;
  
    // k1 = (v1,f1) = F(un) = F(xn,pn) = (vn, q/m E(xn) )
  // v1 = vs = vn
    // calculate_E_mq(a1, xs, xs, q_ws, L, epsilon);
    // v1 = vs;
    if (relativistic) {
        relativistic_momentum_push(xk.data(), v1.data(), pk.data());
    } else {
        nonrelativistic_momentum_push(xk.data(), v1.data(), pk.data());
    }

    for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {

    // k2 = (v2,a2) = F(un + h/2 k1) = F(xn + delt/2 v1, vn + delt/2 a1)
    //              = ( vn + delt a1 / 2, q/m E(xn + delt v1 /2) )
            ef1[xi] *= species_qs[sp_i];

            // k->2
            // p2.push_back(ps[xi] + 0.5 * dt * ef1[xi]);
            pk[xi] = ps[xi] + 0.5 * dt * ef1[xi];
        }
    }
    
    // get x2 = xn + delt v1 / 2
    if (relativistic) {
        relativistic_momentum_push(xk.data(), v2.data(), pk.data());
    } else {
        nonrelativistic_momentum_push(xk.data(), v2.data(), pk.data());
    }
    

    // auto start = high_resolution_clock::now();
    double t2 = t + dt * 0.5;
    evaluate_field(ef2,xk,q_ws,t2);
    // std::vector<double> xk_cpy (xk);
    // std::vector<double> xk_cpy2 (xk);
    // (*calculate_e)(ef2.data(), xk_cpy.data(), ef2.size(),
    //                 xk_cpy2.data(), q_ws.data(), xk.size());
    // auto stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );


    for (size_t sp_i= 0; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            ef2[xi] *= species_qs[sp_i];
            // k -> 3
            pk[xi] = ps[xi] + 0.5 * dt * ef2[xi];
        }
    }
    if (relativistic) {
        relativistic_momentum_push(xk.data(), v3.data(), pk.data());
    } else {
        nonrelativistic_momentum_push(xk.data(), v3.data(), pk.data());
    }


    // k3 = (v3,a3) = F(un + h/2 k2) = F(xn + delt/2 v2, vn + delt/2 a2)
    //              = ( vn + delt a2 / 2, q/m E(xn + delt v2 /2) )
    // for (int ii = 0; ii < N; ++ii) {
    //     xk[ii] = xs[ii] + 0.5 * dt * p2[ii];
    // }

    // start = high_resolution_clock::now();
    double t3 = t2;
    evaluate_field(ef3, xk, q_ws, t3);
    // xk_cpy = xk;
    // xk_cpy2 = xk;
    // (*calculate_e)(ef3.data(), xk_cpy.data(), ef3.size(),
    //                 xk_cpy2.data(), q_ws.data(), xk.size());
    // stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );
    // for (int ii = 0; ii < N; ++ii) {
    //     a3[ii] *= qm;
    // }
    for (size_t sp_i = 0; sp_i < N_sp; ++ sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            ef3[xi] *= species_qms[sp_i];
            // k : 3-> 4
            pk[xi] = ps[xi] + 0.5 * dt * species_qs[sp_i] * ef3[xi];
        }
    }
    // k4 = (v4,a4) = F(un + h k3) = F(xn + delt v3, pn + delt f3)
    //              = ( (pn + delt f3) / m, q E(xn + delt v3) )

    if (relativistic) {
        relativistic_momentum_push(xk.data(), v4.data(), pk.data());
    } else {
        nonrelativistic_momentum_push(xk.data(), v4.data(), pk.data());
    }

    // for (int ii = 0; ii < N; ++ii) {
    //     v4.push_back(vs[ii] + dt * a3[ii]);
    //     xk[ii] = xs[ii] + dt * v3[ii];
    // }
    // start = high_resolution_clock::now();

    double t4 = t + dt;
    evaluate_field(ef4, xk, q_ws, t4);
    // xk_cpy = xk;
    // xk_cpy2 = xk;
    // (*calculate_e)(ef4.data(), xk_cpy.data(), ef4.size(), 
    //                 xk_cpy2.data(), q_ws.data(), xk.size());
    // stop = high_resolution_clock::now();
    // add_time(field_time, duration_cast<duration<double>>(stop - start) );

    // for (int ii = 0; ii < N; ++ii) {
    //     a4[ii] *= qm;
    // }
    for (size_t sp_i = 0; sp_i < N_sp; ++sp_i) {
        for (size_t xi = species_start[sp_i]; xi < species_end[sp_i]; ++xi) {
            ef4[xi] *= species_qms[sp_i];
        }
    }
    // un+1 = (xn+1,pn+1) = un + h/6 (k1 + 2k2 + 2k3 + k4)
    //                    = (xn + delt/6 (v1 + 2 v2 + 2 v3 + v4),
    //                    = pn + delt/6 (a1 + 2 a2 + 2 a3 + a4) )
    // store xn+1, pn+1 in xn,pn
    for (int ii = 0; ii < N; ++ii) {
        xs[ii] += dt / 6.0 * (v1[ii] + 2 * v2[ii] + 2 * v3[ii] + v4[ii]);
        ps[ii] += dt / 6.0 * (ef1[ii] + 2 * ef2[ii] + 2 * ef3[ii] + ef4[ii]);
    }
    if (get_4th_e) {
        // start = high_resolution_clock::now();
        double tk = t + dt;
        evaluate_field(es, xs, q_ws, tk);
        // xk_cpy = xs;
        // xk_cpy2 = xs;
        // (*calculate_e)(es.data(), xk_cpy.data(), xs.size(),
        //                 xk_cpy2.data(), q_ws.data(), xs.size());
        // stop = high_resolution_clock::now();
        // add_time(field_time, duration_cast<duration<double>>(stop - start) );
    }
    

    need_scatter = true;

    return 0;
}
