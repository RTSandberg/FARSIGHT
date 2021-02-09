#include "initial_distributions.hpp"

// #define DEBUG
#ifdef DEBUG
#include <iostream>             // std::cout, std::endl;
using std::cout;
using std::endl;
#include <iterator>             // std::ostream_iterator for vector printing
#endif

// F0_M functions
F0_M::F0_M() : vth(1.0) {}
F0_M::F0_M(double vth) : vth(vth){}

double F0_M::operator() (double x, double v) {
    return 1.0 / sqrt(2.0 * M_PI) / vth * exp(-v*v / 2 / vth /vth);
}
double F0_M::get_vth() { return vth; }
// end F0_M functions

// F0_LD functions
F0_LD::F0_LD() : vth(1.0), k(.5), amp(.1) {}
F0_LD::F0_LD(double vth, double k, double amp) 
    : vth(vth), k(k), amp(amp) {}

double F0_LD::operator() (double x, double v) {
    return 1.0 / sqrt(2.0 * M_PI) / vth * exp(-v*v / 2 / vth /vth) 
        * ( 1 + amp * cos(k * x ));
}
double F0_LD::get_vth() { return vth; }
double F0_LD::get_k() { return k; }
double F0_LD::get_amp() { return amp; }
// end F0_LD functions


// F0_strong_two_stream  functions
F0_strong_two_stream::F0_strong_two_stream() : vth(1.0), k(.5), amp(.1) 
{
    ld_seed = F0_LD(vth, k, amp);
}
F0_strong_two_stream::F0_strong_two_stream(double vth, double k, double amp) 
    : vth(vth), k(k), amp(amp) 
{
    ld_seed = F0_LD(vth, k, amp);
}

double F0_strong_two_stream::operator() (double x, double v) {
    return ld_seed(x,v) * v * v;
}
double F0_strong_two_stream::get_vth() { return vth; }
double F0_strong_two_stream::get_k() { return k; }
double F0_strong_two_stream::get_amp() { return amp; }
// end F0_strong_two_stream functions

// F0_colder_two_stream   functions
F0_colder_two_stream::F0_colder_two_stream() : vth(1.0), v_str(0.5), k(.5), amp(.1) 
{
    maxwellian = F0_M(vth);
}
F0_colder_two_stream::F0_colder_two_stream(double vth, double v_str, double k, double amp) 
    : vth(vth), v_str(v_str), k(k), amp(amp) 
{
    maxwellian = F0_M(vth);
}

double F0_colder_two_stream::operator() (double x, double v) {
    return .5 * (maxwellian(x, v - v_str) 
        + maxwellian(x, v + v_str)) 
        * (1 + amp * cos(k * x));
}
double F0_colder_two_stream::get_vth() { return vth; }
double F0_colder_two_stream::get_k() { return k; }
double F0_colder_two_stream::get_amp() { return amp; }
// end F0_colder_two_stream:: functions

//Friedman beam problem functions
F0_Friedman_beam::F0_Friedman_beam() : Delta(5.522e-8), Tstar(3.463e-7), xmax(0.001)
{
    double vth = sqrt(Tstar);
    nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    mu = 1.25;
    // mu=1.0;
    maxwellian = F0_M(vth);
    generate_Ns(xmax);
}
F0_Friedman_beam::F0_Friedman_beam(double Delta, double Tstar, double xmax) :
    Delta(Delta), Tstar(Tstar), xmax(xmax)
{
    double vth = sqrt(Tstar);
    nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    mu = 1.25;
    // mu=1.0;
    maxwellian = F0_M(vth);
    generate_Ns(xmax);
}
double F0_Friedman_beam::operator() (double x, double v) {
    // return maxwellian(x,v) * nhat * exp(-x*x/2/Lx/Lx);
    return maxwellian(x,v*mu) * nhat * interpolate_N(x/mu);
}
/*
void F0_Friedman_beam::generate_Ns(double xmax) {
    // double xmax = 1.0;
    // convert xmax from kbeta0 xmax to xmax / gammab lambdaD
    double normal_xmax = xmax * x_convert_factor;
    // RK4 integration
    double dx = 0.2;
    int nx = static_cast<int> ( ceil(normal_xmax / dx));
    xs = std::vector<double>(nx);
    Ns = std::vector<double>(nx);
    Ns[0] = 1.0;
    double Np = 0.0;
    for (int ii = 0; ii < nx; ++ii) {
        xs[ii] = ii * dx;
        double Ni = Ns[ii];
        double tN;
        double Np1, Np2, Np3, Np4;
        double dNp1, dNp2, dNp3, dNp4;
        
        Np1 = Np;
        dNp1 = Ni*Ni - Ni * (1 + Delta) + 1/Ni * Np * Np;

        tN = Ni + 0.5 * dx * Np1;
        Np2 = Np + 0.5 * dx * dNp1;
        dNp2 = (tN*tN - tN * (1 + Delta) + 1/tN * Np2* Np2);

        tN = Ni + 0.5 * dx * Np2;
        Np3 = Np + 0.5 * dx * dNp2;
        dNp3 = (tN*tN - tN * (1 + Delta) + 1/tN * Np3 * Np3);

        tN = Ni + 0.5 * dx * Np3;
        Np4 = Np + 0.5 * dx * dNp3;
        dNp4 = (tN*tN - tN * (1 + Delta) + 1/tN * Np4 * Np4);
        
        Ns[ii+1] = Ni + 1.0/6.0 * dx * (Np1 + Np2 + Np3 + Np4);
        Np += 1.0/6.0 * dx * (dNp1 + dNp2 + dNp3 + dNp4);
    }


    // Euler integration
    // double dx = 0.01;
    // int nx = static_cast<int> ( ceil(normal_xmax / dx));
    // xs = std::vector<double>(nx);
    // Ns = std::vector<double>(nx);
    // double Np = 0.0;
    // Ns[0] = 1.0;
    // for (int ii = 0; ii < nx; ++ii) {
    //     xs[ii] = ii * dx;
    //     double Ni = Ns[ii];
    //     Ns[ii+1] = Ni + dx * Np;
    //     Np = Np + dx * (Ni*Ni - Ni * (1 + Delta) + 1/Ni * Np*);
    // }
}
*/
void F0_Friedman_beam::generate_Ns(double xmax) {
    double dx = 0.01;
    double normal_xmax = xmax * x_convert_factor;
    double beyond_xmax = 1.1 * normal_xmax;
    int nx = static_cast<int> ( ceil(beyond_xmax / dx));
    xs = std::vector<double>(nx);
    Ns = std::vector<double>(nx);
    Ns[0] = 1.0;
    double Np = 0.0;
    for (int ii = 1; ii < nx; ++ii) {
        xs[ii] = ii * dx;
        double Nim1 = Ns[ii-1];
        if (Nim1 <= FRIEDMAN_BEAM_MIN_N) { 
            Ns[ii] = 0.0;
        }
        else {
            Ns[ii] = Nim1 + dx * Np;
            Np += dx * (Nim1*Nim1 - Nim1 * (1+Delta) + 1.0/Nim1 * Np*Np);
        }
    }
#ifdef DEBUG
cout << "xs for Friedman beam ic interpolation" << endl;
std::copy(xs.begin(), xs.end(), std::ostream_iterator<double>(cout, " "));
cout << endl;
cout << "Ns for Friedman beam ic interpolation" << endl;
std::copy(Ns.begin(), Ns.end(), std::ostream_iterator<double>(cout, " "));
cout << endl;
#endif /* DEBUG */
}

double F0_Friedman_beam::interpolate_N(double x) {
    double dx = xs[1] - xs[0];
    // convert xmax from kbeta0 xmax to xmax / gammab lambdaD / dx
    double normal_x = fabs(x) * x_convert_factor / dx;
    double ii = floor(normal_x);
    int ind_ii = static_cast<int>(ii);
    double theta = normal_x - ii;
    return (1-theta) * Ns[ind_ii] + theta * Ns[ind_ii+1];
}

//end Friedman beam problem functions
