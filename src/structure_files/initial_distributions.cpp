#include "initial_distributions.hpp"

// #define DEBUG
#ifdef DEBUG
#include <iostream>             // std::cout, std::endl;
using std::cout;
using std::endl;
#include <iterator>             // std::ostream_iterator for vector printing
#endif

// F0_M functions
F0_M::F0_M() : pth(1.0), pstr(0.0) {}
F0_M::F0_M(double pth) : pth(pth), pstr(0.0) {}
F0_M::F0_M(double pth, double pstr) : pth(pth), pstr(pstr) {}

double F0_M::operator() (double x, double p) {
    return 1.0 / sqrt(2.0 * M_PI) / pth * exp(-(p-pstr)*(p-pstr) / 2 / pth /pth);
}
double F0_M::get_pth() { return pth; }
void F0_M::print() {
    cout << "Maxwellian ics" << endl;
    cout << "pth = " << pth << ", pstr = " << pstr << endl;
}
// end F0_M functions

// F0_LD functions
F0_LD::F0_LD() : pth(1.0), pstr(0.0), k(.5), amp(.1) {}
F0_LD::F0_LD(double pth, double k, double amp) 
    : pth(pth), pstr(0.0), k(k), amp(amp) {}
F0_LD::F0_LD(double pth, double pstr, double k, double amp) 
    : pth(pth), pstr(pstr), k(k), amp(amp) {}

double F0_LD::operator() (double x, double p) {
    return 1.0 / sqrt(2.0 * M_PI) / pth * exp(-(p-pstr)*(p-pstr) / 2 / pth /pth) 
        * ( 1 + amp * cos(k * x ));
}
double F0_LD::get_pth() { return pth; }
double F0_LD::get_k() { return k; }
double F0_LD::get_amp() { return amp; }
void F0_LD::print() {
    cout << "Landau damping (sinusoidal perturbation on Maxwellian) ics" << endl;
    cout << "pth = " << pth << ", pstr = " << pstr << ", amp = " << amp << ", k = " << k << endl;
}
// end F0_LD functions


// F0_strong_two_stream  functions
F0_strong_two_stream::F0_strong_two_stream() : pth(1.0), k(.5), amp(.1) 
{
    ld_seed = F0_LD(pth, k, amp);
}
F0_strong_two_stream::F0_strong_two_stream(double pth, double k, double amp) 
    : pth(pth), k(k), amp(amp) 
{
    ld_seed = F0_LD(pth, k, amp);
}

double F0_strong_two_stream::operator() (double x, double p) {
    return ld_seed(x,p) * p * p;
}
double F0_strong_two_stream::get_pth() { return pth; }
double F0_strong_two_stream::get_k() { return k; }
double F0_strong_two_stream::get_amp() { return amp; }
void F0_strong_two_stream::print() {
    cout << "Strong Two stream instability ics" << endl;
    cout << "p^2 * exp (-p^2/2/pth^2) * (1 + amp * cos(k*x))" << endl;
    cout << "pth = " << pth << ", amp = " << amp << ", k = " << k << endl;
}
// end F0_strong_two_stream functions

// F0_colder_two_stream   functions
F0_colder_two_stream::F0_colder_two_stream() : pth(1.0), p_str(0.5), k(.5), amp(.1) 
{
    maxwellian = F0_M(pth);
}
F0_colder_two_stream::F0_colder_two_stream(double pth, double p_str, double k, double amp) 
    : pth(pth), p_str(p_str), k(k), amp(amp) 
{
    maxwellian = F0_M(pth);
}

double F0_colder_two_stream::operator() (double x, double p) {
    return .5 * (maxwellian(x, p - p_str) 
        + maxwellian(x, p + p_str)) 
        * (1 + amp * cos(k * x));
}
double F0_colder_two_stream::get_pth() { return pth; }
double F0_colder_two_stream::get_k() { return k; }
double F0_colder_two_stream::get_amp() { return amp; }
void F0_colder_two_stream::print() {
    cout << "'Colder' Two stream instability ics" << endl;
    cout << "1/2 (fM(p+pstr) + fM(p-pstr)) * (1 + amp * cos(k*x))" << endl;
    cout << "pth = " << pth << ", pstr = " << p_str << ", amp = " << amp << ", k = " << k << endl;
}
// end F0_colder_two_stream:: functions

//Friedman beam problem functions
F0_Friedman_beam::F0_Friedman_beam() : Delta(5.522e-8), Tstar(3.463e-7), xmax(0.001)
{
    double pth = sqrt(Tstar);
    nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    mu = 1.25;
    // mu=1.0;
    maxwellian = F0_M(pth);
    generate_Ns(xmax);
}
F0_Friedman_beam::F0_Friedman_beam(double Delta, double Tstar, double xmax) :
    Delta(Delta), Tstar(Tstar), xmax(xmax)
{
    double pth = sqrt(Tstar);
    nhat = 1.0 / (1.0 + Delta);
    k_beta_0 = 2.0*M_PI / 3.0;
    lambda_D = 2.81e-4;
    beta_b = 2.343e-3;
    gamma_b = 1.0 / sqrt(1 + beta_b * beta_b);
    x_convert_factor = 1/k_beta_0 / lambda_D / gamma_b;
    mu = 1.25;
    // mu=1.0;
    maxwellian = F0_M(pth);
    generate_Ns(xmax);
}
void F0_Friedman_beam::print() {
    cout << "Halo formation ics" << endl;
    cout << "Delta = " << Delta << ", Tstar = " << Tstar << endl;
    cout << "nhat = " << nhat << ", k_beta_0 = " << k_beta_0 << ", lambda_D = " << lambda_D << endl;
    cout << "gamma_b = " << gamma_b << ", mu = " << mu << endl;
}
double F0_Friedman_beam::operator() (double x, double p) {
    // return maxwellian(x,p) * nhat * exp(-x*x/2/Lx/Lx);
    return maxwellian(x,p*mu) * nhat * interpolate_N(x/mu);
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

// Maxwell-Jutner functions
Maxwell_Jutner::Maxwell_Jutner() : theta(1.0), mu(1.0), p_str(0.0), 
                                    k(1.0), amp(0.0), phase_offset(0.0) {}
Maxwell_Jutner::Maxwell_Jutner(double k, double amp) : theta(1.0), mu(1.0), p_str(0.0), 
                                    k(k), amp(amp), phase_offset(0.0) {}
Maxwell_Jutner::Maxwell_Jutner(double mu, double k, double amp) : theta(1.0 / mu), 
                                    mu(mu), p_str(0.0), 
                                    k(k), amp(amp), phase_offset(0.0) {}
Maxwell_Jutner::Maxwell_Jutner(double mu, double p_str, double k, double amp) : 
                                    theta(1.0 / mu), mu(mu), p_str(p_str), 
                                    k(k), amp(amp), phase_offset(0.0) {}
Maxwell_Jutner::Maxwell_Jutner(double mu, double p_str, double k, double amp, double phase_offset) : 
                                    theta(1.0 / mu), mu(mu), p_str(p_str), 
                                    k(k), amp(amp), phase_offset(phase_offset) {}
double Maxwell_Jutner::operator() (double x, double p) {
    double gamma = sqrt(1 + p * p);
    int BesselOrder=1;
    double exp_fac = exp(-gamma * mu) / 2.0 / boost::math::cyl_bessel_k(BesselOrder, mu);
    return exp_fac * (1 + amp * cos(k*x));
    // return 1.0;
}
void Maxwell_Jutner::print() {
    cout << "Perturbation to Relativistic equilibrium ics" << endl;
    cout << "Maxwell-Jutner distribution" << endl;
    cout << "theta = " << theta << ", mu = " << mu << endl;
    cout << "p_str = " << p_str << endl;
    cout << "k = " << k << ", amp = " << amp << ", phase offset = " << phase_offset << endl;
}
// end Maxwell-Jutner functions

// Relativistic two-stream functions
Relativistic_Two_Stream::Relativistic_Two_Stream() 
    : p_str(1.0), mu(1000.0 / 3.0), k(1.0), amp(0.01), phase_offset(0.0) 
{ 
    gamma_str = sqrt(1 + p_str * p_str);
}
Relativistic_Two_Stream::Relativistic_Two_Stream(double p_str, double mu)
    : p_str(p_str), mu(mu), k(1.0), amp(0.01), phase_offset(0.0)
{ 
    gamma_str = sqrt(1 + p_str * p_str);
}
Relativistic_Two_Stream::Relativistic_Two_Stream(double p_str, double mu, double amp)
    : p_str(p_str), mu(mu), k(1.0), amp(amp), phase_offset(0.0)
{ 
    gamma_str = sqrt(1 + p_str * p_str);
}
Relativistic_Two_Stream::Relativistic_Two_Stream(double p_str, double mu, double k, double amp)
    : p_str(p_str), mu(mu), k(k), amp(amp), phase_offset(0.0)
{ 
    gamma_str = sqrt(1 + p_str * p_str);
}
Relativistic_Two_Stream::Relativistic_Two_Stream(double p_str, double mu, double k, double amp, double phase_offset)
    : p_str(p_str), mu(mu), k(k), amp(amp), phase_offset(phase_offset)
{ 
    gamma_str = sqrt(1 + p_str * p_str);
}

double Relativistic_Two_Stream::operator() (double x, double p) {
    double gamma = sqrt(1 + p * p);
    double ppbmu = p * p_str * mu;
    int BesselOrder=1;
    double exp_plus= exp(ppbmu - gamma_str*gamma*mu)
            / 2.0 / boost::math::cyl_bessel_k(BesselOrder, mu);
    double exp_minus = exp(-ppbmu - gamma_str * gamma * mu)
            / 2.0 / boost::math::cyl_bessel_k(BesselOrder, mu);
    return (exp_plus + exp_minus) / 2.0 * (1 + amp*cos(k*x));
}
void Relativistic_Two_Stream::print()  {
    cout << "Relativistic two-stream instability initial conditions" << endl;
    cout << "Two relativistic Maxwellians" << endl;
    cout << "mu = " << mu << ", p_str = " << p_str << endl;
    cout << "k = " << k << ", amp = " << amp << ", phase offset = " << phase_offset << endl;
}
// end Relativistic two-stream functions

// Relativistic wave functions
int check_gamma(double max_gamma, double wave_beta) {

    if (wave_beta < 1.0) {
        double gamma_wave = 1.0 / sqrt(1 - wave_beta*wave_beta);
        if (gamma_wave < max_gamma) {
            throw std::invalid_argument("Can't have particles moving faster than wave! Choose smaller max gamma or faster wave");
        }
    }
    return 0;
}

Relativistic_Wave::Relativistic_Wave() 
    : amp(1.01), wave_beta(0.8), p_th(1.0) 
{   
    make_cold_xs_ps(); 
}
Relativistic_Wave::Relativistic_Wave(double amp)
    : amp(amp), wave_beta(0.8), p_th(1.0) 
{ 
    check_gamma(amp, wave_beta);
    make_cold_xs_ps(); 
}
Relativistic_Wave::Relativistic_Wave(double amp, double wave_beta)
    : amp(amp), wave_beta(wave_beta), p_th(1.0) 
{   
    check_gamma(amp, wave_beta);
    make_cold_xs_ps(); 
}
Relativistic_Wave::Relativistic_Wave(double amp, double wave_beta, double p_th)
    : amp(amp), wave_beta(wave_beta), p_th(p_th) 
{   
    check_gamma(amp, wave_beta);
    make_cold_xs_ps(); 
}
void Relativistic_Wave::make_cold_xs_ps() {
    int Np = 100;
    double gm = amp;
    double kappa_sq = (gm-1.0) / (gm + 1.0);
    double kappa = sqrt(kappa_sq);
    double kappa_p = sqrt(2 / (gm + 1.0));

    xs = std::vector<double> (Np);
    ps = std::vector<double> (Np);
    double dalpha = 3*M_PI / (Np - 1);
    try {
    for (int ii = 0; ii < Np; ++ii) { 
        double alphai = -M_PI/2 + ii * dalpha;
        double x0i = wave_beta * 
                (2 / kappa_p * boost::math::ellint_2(kappa, alphai)
                    - kappa_p * boost::math::ellint_1(kappa, alphai));
        double dxi = 2 * kappa / kappa_p * sin(alphai);
        xs[ii] = x0i + dxi;

        double gammai = gm - (gm - 1) * sin(alphai)*sin(alphai);
        double abs_pi = sqrt(gammai*gammai - 1);
        double signp = (cos(alphai)) > 0 ? 1 : -1;
        ps[ii] = abs_pi * signp;
    }
    }
    catch(const std::exception& e) {
        cout <<
         "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
    }
}
double Relativistic_Wave::operator() (double x, double p) {
    // interpolate from cold wave solution
    int ind_ii = 0;
    while (x > xs[ind_ii+1]) {
        ++ind_ii;
    }
    double theta = x - xs[ind_ii];
    double p0 = (1-theta) * ps[ind_ii] + theta * ps[ind_ii+1];
    
    // if (fabs(p-p0) < p_th) {
    //     return 1.0/2.0/p_th;
    // }
    // else {
    //     return 0;
    // }

    return 1.0/sqrt(2*M_PI)/p_th*exp(-(p-p0)*(p-p0)/2/p_th/p_th);
}
void Relativistic_Wave::print()  {
    cout << "Relativistic travelign wave initial conditions" << endl;
    cout << "Thermal width applied to cold solution" << endl;
    cout << "p_th = " << p_th << endl;
    cout << "wave beta = " << wave_beta << ", amp = " << amp <<  endl;
}
// end relativistic wave functions
