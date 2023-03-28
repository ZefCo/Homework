#include <cmath>
#include <vector>
#include <array>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <functional>
#include <tuple>

namespace fs = std::filesystem;

const double hbar = 1.0546E-34;  // hbar
const double me = 9.109E-31;     // mass of electron
const double e = 1.6022E-19;     // electron charge
const double a = 10E-11;         // Bohr radius (order of magnitude)



class ZeroWell {
    public:
        // constructor
        ZeroWell(double V0, double Ll, double Lr, double N, double E1, double E2);
        // destructor
        ~ZeroWell();

        // uses RK4 to find the boundary values
        std::tuple<double, double> wavefunction(double e);

        // it may seem silly but this has the same form as the other two so I can test the rest of the code
        // and be lazy at the same time: don't have to alter the inputs to the function.
        // Ideally I would instead pass in the potential as a class parameter, but I'm not competant enough 
        // in C++ to do that
        double V(double x, double a);

        std::tuple<double, double> func(double psi, double phi, double x, double E);

        void secant_ODE(double target);

        double get_E2();




    private:
        // V0 of the potential
        double V0;
        // the left and right end points of the well
        double Ll, Lr;
        // number of points in the well to calculate at
        double N;
        // step size of the well
        double h;
        // Energy of the particle
        double E1;
        // Guess energy of the particle
        double E2;

        // the vector that stores the points to caluclate the wavefunction at
        std::vector<double> well_points;

        // Creates a vector of all the points for the wavefunction to be calculated at. I don't want to redo
        // this over and over again like in the example code. It puts these values in the well_points vector.
        void potential_points();

        // double potential(double, double, double);

};


class HarmonicWell {
    public:
        // constructor
        HarmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2);
        // destructor
        ~HarmonicWell();

        // uses RK4 to find the boundary values
        std::tuple<double, double> wavefunction(double e);

        // A potential of V(x) = V0 x**2 / a**2
        double V(double x, double a);

        std::tuple<double, double> func(double psi, double phi, double x, double E);

        void secant_ODE(double target);

        double get_E2();




    private:
        // V0 of the potential
        double V0;
        // the left and right end points of the well
        double Ll, Lr;
        // number of points in the well to calculate at
        double N;
        // step size of the well
        double h;
        // Energy of the particle
        double E1;
        // Guess energy of the particle
        double E2;

        // the vector that stores the points to caluclate the wavefunction at
        std::vector<double> well_points;

        // Creates a vector of all the points for the wavefunction to be calculated at. I don't want to redo
        // this over and over again like in the example code. It puts these values in the well_points vector.
        void potential_points();

        // double potential(double, double, double);

};



class AnharmonicWell {
    public:
        // constructor
        AnharmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2);
        // destructor
        ~AnharmonicWell();

        // uses RK4 to find the boundary values
        std::tuple<double, double> wavefunction(double e);

        // A potential of V(x) = V0 x**4 / a**4
        double V(double x, double a);

        std::tuple<double, double> func(double psi, double phi, double x, double E);

        void secant_ODE(double target);

        double get_E2();




    private:
        // V0 of the potential
        double V0;
        // the left and right end points of the well
        double Ll, Lr;
        // number of points in the well to calculate at
        double N;
        // step size of the well
        double h;
        // Energy of the particle
        double E1;
        // Guess energy of the particle
        double E2;

        // the vector that stores the points to caluclate the wavefunction at
        std::vector<double> well_points;

        // Creates a vector of all the points for the wavefunction to be calculated at. I don't want to redo
        // this over and over again like in the example code. It puts these values in the well_points vector.
        void potential_points();

        // double potential(double, double, double);

};