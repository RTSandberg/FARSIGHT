#include "FieldStructure.hpp"

ElectricField::~ElectricField() = default;


E_MQ_DirectSum::E_MQ_DirectSum() {}
E_MQ_DirectSum::E_MQ_DirectSum(double L, double epsilon) : L(L), epsilon(epsilon) {}
E_MQ_DirectSum::~E_MQ_DirectSum() = default;

void E_MQ_DirectSum::operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
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
}

E_MQ_DirectSum_openbcs::E_MQ_DirectSum_openbcs() {}
E_MQ_DirectSum_openbcs::E_MQ_DirectSum_openbcs(double epsilon) : epsilon(epsilon) {}
E_MQ_DirectSum_openbcs::~E_MQ_DirectSum_openbcs() = default;

void E_MQ_DirectSum_openbcs::operator() (double* es, double* targets, int nt, 
                        double* sources, double* q_ws, int ns)
{    
    double eps_sq = epsilon * epsilon;

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
            double z = xi - sources[jj];
            ei += q_ws[jj] * 0.5 * z / sqrt(z * z + eps_sq);
        }
        es[ii] = ei - xi;
    }
}

E_MQ_Treecode::E_MQ_Treecode() {}
E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon, double beta) : 
    kernel(MQ), 
    singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
    beta(beta), theta(-1.0), interpDegree(-1), maxPerSourceLeaf(1), maxPerTargetLeaf(1),
    verbosity(0)
{
    kernelParams.push_back(L); kernelParams.push_back(epsilon);
}

E_MQ_Treecode::E_MQ_Treecode(double L, double epsilon,
    double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
    int verbosity) :
        kernel(MQ), 
        singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
        beta(-1.0), theta(theta), interpDegree(interpDegree), 
        maxPerSourceLeaf(maxPerSourceLeaf), maxPerTargetLeaf(maxPerTargetLeaf),
        verbosity(verbosity)
{
    kernelParams.push_back(L); kernelParams.push_back(epsilon);
}

void E_MQ_Treecode::operator() (double* es, double* targets, int nt, 
                double* sources, double* q_ws, int ns)
{
    std::vector <double> xS(ns);
    std::vector <double> yS(ns);
    std::vector <double> wS(ns, 1.0);

    std::vector <double> xT(nt);
    std::vector <double> yT(nt);
    std::vector <double> qT(nt, 1.0);

    // if (targets == sources) {
    //     std::vector<double> new_sources(ns);
    //     for (int ii = 0; ii < ns; ii++) {
    //         new_sources[ii] = sources[ii];
    //     }
    //     sources = new_sources.data();
    // }
//
/*
    std::cout << "beta " << beta << std::endl;
    std::cout << "mac " << theta << std::endl;
    std::cout << "degree " << interpDegree << std::endl;
    std::cout << "max source " << maxPerSourceLeaf << std::endl;
    std::cout << "max target " << maxPerTargetLeaf << std::endl;
    std::cout << "ns " << ns << std::endl;
    std::cout << "nt " << nt << std::endl;
    std::cout << "kernel " << kernel << std::endl;
    std::cout << "singularity " << singularity << std::endl;
    std::cout << "approximation " << approximation << std::endl;
    std::cout << "compute type " << compute_type << std::endl;

    std::cout << "xT ";
    std::copy(xT.begin(), xT.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "xS ";
    std::copy(xS.begin(), xS.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "yT ";
    std::copy(yT.begin(), yT.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "yS ";
    std::copy(yS.begin(), yS.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "qT ";
    std::copy(qT.begin(), qT.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "wS ";
    std::copy(wS.begin(), wS.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "es ";
    for (int ii = 0; ii < nt; ii++) {
        std::cout << es[ii] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "xs ";
    for (int ii = 0; ii < ns; ii++) {
        std::cout << sources[ii] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "xt ";
    for (int ii = 0; ii < nt; ii++) {
        std::cout << targets[ii] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "qws ";
    for (int ii = 0; ii < ns; ii++) {
        std::cout << q_ws[ii] << " ";
    }
    std::cout << std::endl;   
    */

    BaryTreeInterface(nt, ns, xT.data(), yT.data(), 
                    targets, qT.data(),
                    xS.data(), yS.data(), 
                    sources, q_ws, wS.data(),
                    es,
                    kernel, kernelParams.size(), kernelParams.data(),
                    singularity, approximation, compute_type,
                    theta, interpDegree, maxPerSourceLeaf, maxPerTargetLeaf,
                    1.0, beta, verbosity);
}
    
E_MQ_Treecode::~E_MQ_Treecode() = default;

