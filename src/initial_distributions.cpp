#include "initial_distributions.hpp"

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
F0_Friedman_beam::F0_Friedman_beam() : nhat(0.99999994478), Tstar(3.463e-7), xmax(0.001)
{
    double vth = sqrt(Tstar);
    // double nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    maxwellian = F0_M(vth);
    generate_Ns(xmax);
}
F0_Friedman_beam::F0_Friedman_beam(double nhat, double Tstar, double xmax) :
    nhat(nhat), Tstar(Tstar), xmax(xmax)
{
    double vth = sqrt(Tstar);
    // nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    maxwellian = F0_M(vth);
    generate_Ns(xmax);
}
double F0_Friedman_beam::operator() (double x, double v) {
    // return maxwellian(x,v) * nhat * exp(-x*x/2/Lx/Lx);
    return maxwellian(x,v) * nhat * interpolate_N(x);
}
void F0_Friedman_beam::generate_Ns(double xmax) {
    // double xmax = 1.0;
    // convert xmax from kbeta0 xmax to xmax / gammab lambdaD
    double normal_xmax = xmax * x_convert_factor;
    double dx = 0.01;
    int nx = static_cast<int> ( ceil(normal_xmax / dx));
    xs = std::vector<double>(nx);
    Ns = std::vector<double>(nx);
    double Np = 0.0;
    for (int ii = 0; ii < nx; ++ii) {
        xs[ii] = ii * dx;
        double Ni = Ns[ii];
        Ns[ii+1] = Ni + dx * Np;
        Np = Np + dx * (Ni*Ni - Ni * (1 + Delta) + 1/Ni * Np);
    }
}
double F0_Friedman_beam::interpolate_N(double x) {
    // convert xmax from kbeta0 xmax to xmax / gammab lambdaD
    double normal_x = x * x_convert_factor;
    double dx = xs[1] - xs[0];
    double ii = floor(fabs(x) / dx);
    int ind_ii = static_cast<int>(ii);
    double theta = fabs(normal_x) - ii;
    return (1-theta) * Ns[ii] + theta * Ns[ii+1];
}

//end Friedman beam problem functions
