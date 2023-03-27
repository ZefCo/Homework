#include <vector>
#include <array>
#include <filesystem>
#include <iostream>
#include <fstream>

namespace fs = std::filesystem;

// it may seem silly but this has the same form as the other two so I can test the rest of the code
// and be lazy at the same time: don't have to alter the inputs to the function
double V_zero(double x, double a, double V0);

// A potential of V(x) = V0 x**2 / a**2
double V_harmonic_oscillator(double x, double a, double V0);

// A potential of V(x) = V0 x**4 / a**4
double V_anharmoinc_oscillator(double x, double a, double V0);


// uses RK4 to find the boundary values
void wavefunction();


// Returns a vector of all the points for the wavefunction to be calculated at. I don't want to redo
// this over and over again like in the example code
std::vector<double> potential_points(double start, double stop, double step);