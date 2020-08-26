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
    double beta, int verbosity) :
        kernel(MQ), 
        singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
        beta(beta), theta(theta), interpDegree(interpDegree), 
        maxPerSourceLeaf(maxPerSourceLeaf), maxPerTargetLeaf(maxPerTargetLeaf),
        verbosity(verbosity) 
{
    kernelParams.push_back(L); kernelParams.push_back(epsilon);
}

void E_MQ_Treecode::operator() (double* es, double* targets, int nt, 
                double* sources, double* q_ws, int ns)
{
    std::vector <double> xS (ns);
    std::vector <double> yS (ns);
    std::vector <double> wS (ns);

    std::vector <double> xT (nt);
    std::vector <double> yT (nt);
    std::vector <double> qT (nt);

    std::fill(wS.begin(), wS.end(), 1.);
    std::fill(qT.begin(), qT.end(), 1.);

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

