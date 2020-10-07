#include "AMRStructure.hpp"

int calculate_E_mq(double* es, const double* targets, int nt, 
                const double* sources, const double* q_ws, int ns, 
                double L, double epsilon) {

    double epsLsq = epsilon * epsilon / L / L;
    double norm_epsL = sqrt(1 + 4 * epsLsq );

#ifdef OPENACC_ENABLED
#pragma acc parallel loop independent
#else
#pragma omp parallel for
#endif
    for (int ii = 0; ii < nt; ++ii) {
        double xi = targets[ii];
        double ei = 0.0;
#ifdef OPENACC_ENABLED
#pragma acc loop independent reduction(+:ei)
#endif
        for (int jj = 0; jj < ns; ++jj) {
            double z = (xi - sources[jj]) / L;

            while (z < -0.5) { z += 1.0; }
            while (z >= 0.5) { z -= 1.0; }
            ei += q_ws[jj] * (0.5 * z * norm_epsL / sqrt(z * z + epsLsq) - z);
        }
        es[ii] = ei;
    }
    return 0;
}

// void AMRStructure::vector_e_wrapper() {

//     std::vector<double> xs, es, qws;
//     xs.reserve(particles.size());
//     es.reserve(particles.size());
//     qws.reserve(particles.size());
//     for (auto& particle : particles) {
//         xs.push_back(particle.x);
//         qws.push_back(particle.q_w);
//         es.push_back(particle.e);
//     }
//     calculate_E_mq(es, xs, xs, qws, L, epsilon);
//     for(int ii = 0; ii < es.size(); ii++) {
//         particles[ii].e = es[ii];
//     }
// }
// void AMRStructure::vector_calc_e_wrapper() {
//     // es.reserve(xs.size());
//     es = std::vector<double>(xs.size());
//     calculate_E_mq(es.data(), xs.data(), xs.size(), 
//                 xs.data(), q_ws.data(), xs.size(),
//                 Lx, greens_epsilon);
// }


void AMRStructure::init_e() {
    // get unique sort xs, inv_inds
    std::vector<size_t> inds(xs.size()); //(arr_inds, arr_inds +  sizeof(arr_inds) / sizeof(size_t));
    std::iota(inds.begin(), inds.end(), 0);
    std::sort(inds.begin(), inds.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });

    // unique pat
    int s_ind = inds[0];
    double val = xs[s_ind];
    std::vector<double> unique_xs (1,val);
    std::vector<size_t> inv_inds (xs.size());
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
            unique_xs.push_back(xsi);
        } 
        // std::cout << "inv_inds[" << s_ind << "] set to " << unique_ind << std::endl;
        inv_inds[s_ind] = unique_ind;
    }

    std::vector<double> sort_ws (unique_xs.size());
    for (int ii = 0; ii < inv_inds.size(); ++ii) {
        sort_ws[inv_inds[ii]] += q_ws[ii];
    }
    // calculate reduced E
    std::vector<double> sort_es(sort_ws.size());

    // #pragma omp parallel 
    // {
    // int num_threads = omp_get_num_threads();
    // if (omp_get_thread_num() == 0) {
        // cout << "Number of threads from init_e(): " << num_threads << endl;
    // }
    // calculate_E_mq(sort_es.data(), unique_xs.data(), unique_xs.size(),
    //                 unique_xs.data(), sort_ws.data(), unique_xs.size(),
    //                 Lx, greens_epsilon);
    // calculate_e(sort_es.data(), )


    std::vector<double> unique_xs_cpy (unique_xs);

    (*calculate_e)(sort_es.data(), unique_xs.data(), unique_xs.size(),
                    unique_xs_cpy.data(), sort_ws.data(), unique_xs.size() );
    // }
    bool do_subtract_mean_e = false;
    if (do_subtract_mean_e) {
        double mean_e = 0;
        for (int ii = 0; ii < sort_es.size(); ++ii) {
            mean_e += sort_ws[ii] * sort_es[ii];
        }
        mean_e *= 1.0 / Q0;
        for (int ii = 0; ii < sort_es.size(); ++ii) {
            sort_es[ii] -= mean_e;
        }
    }
    // es.clear();
    es = std::vector<double>(xs.size());
    for (int ii = 0; ii < inv_inds.size(); ++ii) {
        // es.push_back(sort_es[inv_inds[ii]]);
        es[ii] = sort_es[inv_inds[ii]];
    }
}

// void AMRStructure::calculate_e() {
//     double epsLsq = epsilon * epsilon / L / L;
//     double norm_epsL = sqrt(1 + 4 * epsLsq );

//     for (int ii = 0; ii < particles.size(); ++ii) {
//         double xi = particles[ii].x;
//         for (int jj = 0; jj < particles.size(); ++jj) {
//             double z = (xi - particles[jj].x) / L;

//             while (z < -0.5) { z += 1.0; }
//             while (z >= 0.5) { z -= 1.0; }
//             double eij = particles[jj].q_w * (0.5 * z * norm_epsL / sqrt(z * z + epsLsq) - z);
//             particles[ii].e += eij;
//         }
//     }
// }
