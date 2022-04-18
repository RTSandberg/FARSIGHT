#include "FieldStructure.hpp"

ExternalElectricField::~ExternalElectricField() = default;
// -----------------
ExternalTanh::ExternalTanh()
        : Am(0.2), k(0.26), omega(0.37) {
            eps = drive_fun(t0);
        }
ExternalTanh::ExternalTanh(double Am)
        : Am(Am), k(0.26), omega(0.37) {
            eps = drive_fun(t0);
        }    
ExternalTanh::ExternalTanh(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {
            eps = drive_fun(t0);
        }   
double ExternalTanh::drive_fun(double t) {
    return  0.5 * (tanh((t - tL)/twL) - tanh((t-tR)/twR));
}
void ExternalTanh::operator() (double* es, double* targets, int nt, double t) {
    double a_t_k = (drive_fun(t) - eps) / (1-eps) * Am * k;
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a_t_k * sin(k*targets[ii] - omt);
    }
}
void ExternalTanh::print_field_obj() {
    cout << "External field : Tanh" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
    cout << "t0 " << t0 << ", tL " << tL << ", tR " << tR << endl;
    cout << "twL " << twL << ", twR " << twR << endl;
    cout << "epsilon driver " << eps << endl;
}
ExternalTanh::~ExternalTanh() = default;
// -----------------
ExternalSine::ExternalSine()
        : Am(0.4), k(0.26), omega(0.37) {}
ExternalSine::ExternalSine(double Am)
        : Am(Am), k(0.26), omega(0.37) {}    
ExternalSine::ExternalSine(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {}   
void ExternalSine::operator() (double* es, double* targets, int nt, double t) {
    double a0 = 0;
    if (t < t1) {
        a0 = Am * sin (t * M_PI / 100);
    } else if (t < t2) {
        a0 = Am;
    } else if (t < t3) {
        a0 = Am * cos((t-150)*M_PI / 100);
    }
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a0 * sin(k*targets[ii] - omt);
    }
}
void ExternalSine::print_field_obj() {
    cout << "External field : sine" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
}
ExternalSine::~ExternalSine() = default;
//--------------------------
ExternalLogistic::ExternalLogistic()
        : Am(0.4), k(0.26), omega(0.37) {}
ExternalLogistic::ExternalLogistic(double Am)
        : Am(Am), k(0.26), omega(0.37) {}    
ExternalLogistic::ExternalLogistic(double Am, double k, double omega)
        : Am(Am), k(k), omega(omega) {}   
void ExternalLogistic::operator() (double* es, double* targets, int nt, double t) {
    double a0;
    if (t < t1) {
        a0 = Am / (1.0 + exp(-40.0*(t-10.0)));
    } else {
        a0 = Am * (1.0 - 1.0 / (1.0 + exp(-40.0 * (t - 110.0))));
    }
    double omt = omega * t;
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] += a0 * sin(k*targets[ii] - omt);
    }
}
void ExternalLogistic::print_field_obj() {
    cout << "External field : logistic" << endl;
    cout << "Am " << Am << ", k " << k << ", omega " << omega << endl;
}
ExternalLogistic::~ExternalLogistic() = default;
//================== end external ======

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
void E_MQ_DirectSum::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, direct sum" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "Periodic boundary conditions, domain size = " << L << endl;
    cout << "-------------" << endl;
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
void E_MQ_DirectSum_openbcs::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, direct sum" << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "Open boundary conditions "<< endl;
    cout << "-------------" << endl;
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
void E_MQ_Treecode::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, treecode" << endl;
    cout << "epsilon = " << kernelParams[1] << endl;
    cout << "Periodic boundary conditions, domain size = " << kernelParams[0] << endl;
    if (beta >= 0) {
        cout << "treecode parameters set with beta = " << beta << endl;
    } else {
        cout << "mac = " << theta << ", n = " << interpDegree << endl;
        cout << " max source = " << maxPerSourceLeaf << " , max target = " << maxPerTargetLeaf << endl;
    }
    cout << "-------------" << endl;
}
    
E_MQ_Treecode::~E_MQ_Treecode() = default;


E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs() {}
E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs(double epsilon, double beta) : 
    kernel(MQ), 
    singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
    beta(beta), theta(-1.0), interpDegree(-1), maxPerSourceLeaf(1), maxPerTargetLeaf(1),
    verbosity(0)
{
    kernelParams.push_back(-1.0); kernelParams.push_back(epsilon);
}

E_MQ_Treecode_openbcs::E_MQ_Treecode_openbcs(double epsilon,
    double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
    int verbosity) :
        kernel(MQ), 
        singularity(SKIPPING), approximation(LAGRANGE), compute_type(PARTICLE_CLUSTER),
        beta(-1.0), theta(theta), interpDegree(interpDegree), 
        maxPerSourceLeaf(maxPerSourceLeaf), maxPerTargetLeaf(maxPerTargetLeaf),
        verbosity(verbosity)
{
    kernelParams.push_back(-1.0); kernelParams.push_back(epsilon);
}

void E_MQ_Treecode_openbcs::operator() (double* es, double* targets, int nt, 
                double* sources, double* q_ws, int ns)
{
    std::vector <double> xS(ns);
    std::vector <double> yS(ns);
    std::vector <double> wS(ns, 1.0);

    std::vector <double> xT(nt);
    std::vector <double> yT(nt);
    std::vector <double> qT(nt, 1.0);

    BaryTreeInterface(nt, ns, xT.data(), yT.data(), 
                    targets, qT.data(),
                    xS.data(), yS.data(), 
                    sources, q_ws, wS.data(),
                    es,
                    kernel, kernelParams.size(), kernelParams.data(),
                    singularity, approximation, compute_type,
                    theta, interpDegree, maxPerSourceLeaf, maxPerTargetLeaf,
                    1.0, beta, verbosity);
    for (int ii = 0; ii < nt; ++ii) {
        es[ii] -= targets[ii];
    }
}
void E_MQ_Treecode_openbcs::print_field_obj() {
    cout << "-------------" << endl;
    cout << "Field object: " << endl;
    cout << "MQ kernel, treecode" << endl;
    cout << "epsilon = " << kernelParams[0] << endl;
    cout << "Open boundary conditions" << endl;
    if (beta >= 0) {
        cout << "treecode parameters set with beta = " << beta << endl;
    } else {
        cout << "mac = " << theta << ", n = " << interpDegree << endl;
        cout << " max source = " << maxPerSourceLeaf << " , max target = " << maxPerTargetLeaf << endl;
    }
    cout << "-------------" << endl;
}
    
E_MQ_Treecode_openbcs::~E_MQ_Treecode_openbcs() = default;
