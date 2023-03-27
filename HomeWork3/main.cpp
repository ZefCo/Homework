// #include <iostream>
#include <tuple>
// #include <filesystem>
// #include "json/json.h"
#include "Problem1.h"
#include "Problem2.h"
#include "Problem3.h"

// Windows
// g++ main.cpp Problem1.cpp Problem2.cpp Problem3.cpp -o HW3.exe
// Linux
// g++ main.cpp Problem1.cpp Problem2.cpp Problem3.cpp -o HW3.out

namespace fs = std::filesystem;

// Imported in Problem 2
// const double hbar = 1.0546E-34;  // hbar
// const double me = 9.109E-32;     // mass of electron
// const double e = 1.6022E-19;     // electron charge
// const double a = 10E-11;         // Bohr radius (order of magnitude)


void problem1() {
    double s = 10, r = 28, b = 8/3, h = 0.01;
    double tfinal = 50;
    double x0 = 0, y0 = 1, z0 = 0, t0 = 0;

    std::cout << "Problem 1" << std::endl;

    fs::path outfile1a = fs::current_path() / "Lorenz_Problem_1a.csv";

    std::cout << "Values for Problem 1 a)" << std::endl;
    std::cout << "x0 = " << x0 << " y0 = " << y0 << " z0 = " << z0 << std::endl;
    std::cout << "s = " << s << " r = " << r << " b = " << b << std::endl;
    std::cout << "t0 = " << t0 << " t final = " << tfinal << " t size = " << h << std::endl;

    LorenzRK4 lorenz_y(x0, y0, z0, t0, tfinal, s, r, b, h);
    lorenz_y.ODE_Solver();

    lorenz_y.write_csv(outfile1a);

}


void problem2() {
    double v0, L, N, h;

    v0 = (50)*e;         // puts V0 into eV
    L = 10*a - (-10*a);  // size of the well
    N = 1000;            // number of points in the well
    h = L / N;           // step size of the well

    // std::vector<double> ppoints = potential_points(-10*a, 10*a, h);


    // std::cout << "L = " << L << " N = " << N << " a = " << a << " h = " << h << std::endl;

    // std::cout << "Size of ppoints = " << ppoints.size() << std::endl;
    // std::cout << "First three points of ppoint:\n" << ppoints[0] << ", " << ppoints[1] << ", " << ppoints[2] << std::endl;
    // std::cout << "Last three points of ppoint:\n" << ppoints[ppoints.size() - 3] << ", " << ppoints[ppoints.size() - 2] << ", " << ppoints[ppoints.size() - 1] << std::endl;
    
    ZeroWell test(v0, -10*a, 10*a, N, 0, e);

    test.secant_ODE(e/1000);

    std::cout << "Part a E2 = " << test.get_E2() / e << std::endl;

}


int main() {
    // std::cout << "hbar = " << hbar << std::endl;

    std::cout << "Problem 1 turned off right now" << std::endl;
    // problem1();
    problem2();


}