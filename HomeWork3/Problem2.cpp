#include "Problem2.h"

// ParticleInABox::ParticleInABox()


double V_zero(double x, double a, double V0) {
    return 0.0;
}

double V_harmonic_oscillator(double x, double a, double V0) {
    return V0 * x*x / (a*a);
}

double V_anharmoinc_oscillator(double x, double a, double V0) {
    return V0 * (x*x*x*x) / (a*a*a*a);
}

void wavefunction() {

}


std::vector<double> potential_points(double start, double stop, double step) {
    std::vector<double> return_vec;

    for (double i = start; i = stop; i+= step) {
        return_vec.push_back(i);
    }

    return return_vec;

}