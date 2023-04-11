#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;


// I probably should fold these into the class... but I wont
const double tau = 365;                       // years
const double A = 10;                          // degrees C
const double B = 12;                          // degrees C
const double D = 0.1;                         // m**2 day**-1
const double PI = 3.141592653589793238463;    // value of pi, which I'm kind of surprised is not built into cmath or math.h...


class EarthCrust {
    public:
        EarthCrust(double years, double TEarthCore, double TEarthSurface, int N, int L, double h);
        ~EarthCrust();

        // function for the tempurature
        double T(double t);

        // stores the values for the FTCS method. I'm just adapting the one in the book, I'm not going to rethink and try to improve upon the algorithm.
        // I've got an APE in like a month...
        std::vector<double> depth;
        std::vector<double> depth_prime;

        // Runs the simulation for N years. Runs it for N - 1 Years, then outputs the data for every three months (1 quarter) in the Nth year.
        // So if you did 10 years it would run for 9, then output the last year in 4 chunks.
        std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> run_sim();

    private:
        // Yes I know the Earth's core is much, much more then 11 degrees C, but it's just easier to call it this
        double TEarthCore;

        // Surface value, also used as the intermediate value
        double TEarthSurface;

        // Number of points on the grid
        int N;
        
        // thickness of the Earth/length of interest
        int L;

        // Creates a vector with index 0 as the TEarthCore value, and TEarthSurface for all other values
        void init_grid();

        // To save on compute time: c = h * D / a**2, the constant sitting outside the FTCS (9.20 from book)
        double c;

        // the step size
        double h;

        // the grid spacing
        double a;

        // a buffer value that goes a small step beyond, which allows me to catch each quarter as they process.
        double epsilon;

        // final time point for the simulation. User imputs the years, this figures it out in days.
        double tend;

        // it's easier to use a double then to use an int here
        double years;

        // values for when to switch to the different sections. Warmup isn't really used, but the rest are part of an if statement: 
        // if |t - t'|< epsilon then the vector that was just computed represents the given quarter of interest and that is saved to a special vector that is output from the run simulation.
        double warmup;
        double first_quarter;
        double second_quarter;
        double third_quarter;
        double fourth_quarter;

};


// Writes the four vectors that come from run_sim() to a csv file.
void write_problem1(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d);