// #include "json/json.h"
#include <vector>
#include <array>
#include <filesystem>
#include <iostream>
#include <fstream>

namespace fs = std::filesystem;

const int X = 0, Y = 1, Z = 2;


class LorenzRK4 {
    public:
        // constructor
        LorenzRK4(double x0, double y0, double z0, double t0, double tfinal, double sigma, double r, double b, double h);
        // destructor
        ~LorenzRK4();

        // Uses the RK4 method to solve the equations
        void ODE_Solver();

        // does a single step of the RK4 method then pushes the results to the xn, yn, & zn vectors
        void time_step(int t);

        // dx/dt = s(y - x)
        double dxdt(double x, double y, double z);
        // dy/dt = r*x - y- x*z
        double dydt(double x, double y, double z);
        // dz/dt = x*y - b*z
        double dzdt(double x, double y, double z);
        // keeps track of the time
        void time(int t); 

        // creates a blank table
        void init_table();

        // writes the results to a csv file
        void write_csv(fs::path outfile);

    private:

        double s;
        double r;
        double b;
        double h;
        double tfinal;
        int tmax;

        std::vector<double> xn;
        std::vector<double> yn;
        std::vector<double> zn;
        std::vector<double> tn;

        double x_store;
        double y_store;
        double z_store;

        // an empty table: I can just call this intstead of recreating a blank table over and over again
        std::vector<std::vector<double>> blank_table;

};
