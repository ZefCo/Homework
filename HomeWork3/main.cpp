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


void problem1() {
    const double s = 10, r = 28, b = 8/3, h = 0.01;
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

    std::cout << "Problem 1 printed to csv: use HW3_plots.ipynb for graphs" << std::endl << std::endl;

}


void problem2() {
    double v0, L, N, h;

    std::cout << "Problem 2" << std::endl;

    v0 = (50)*e;         // puts V0 into eV
    N = 1000;            // number of points in the well
    h = L / N;           // step size of the well
    // other constants defined in Problem2.h file
    
    ZeroWell test(v0, 0, 5.2918E-11, N, 0, e);
    test.secant_ODE(e/1000);
    std::cout << " E/e = " << test.get_E2()/e << "ev" << std::endl;


    std::cout << "Part a)" << std::endl;
    HarmonicWell partb_g(v0, -10*a, 10*a, N, 0, e);
    partb_g.secant_ODE(e/1000.0);
    double partb_ge = partb_g.get_E2() / e;
    std::cout << "E g state = " << partb_ge << " ev" << std::endl;

    HarmonicWell partb_1(v0, -10*a, 10*a, N, partb_ge*e, 3*partb_ge*e);
    partb_1.secant_ODE(e/1000.0);
    double partb_1e = partb_1.get_E2() / e;
    std::cout << "E 1 state = " << partb_1e << " ev" << std::endl;

    HarmonicWell partb_2(v0, -10*a, 10*a, N, 4*partb_ge*e, 7*partb_ge*e);
    partb_2.secant_ODE(e/1000.0);
    double partb_2e = partb_2.get_E2() / e;
    std::cout << "E 2 state = " << partb_2e << " ev" << std::endl;


}


void problem3() {

}


int main() {
    // std::cout << "hbar = " << hbar << std::endl;

    std::cout << "Problem 1 turned off right now" << std::endl;
    // problem1();
    problem2();


}