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

    std::cout << "Problem 1\n" << std::endl;

    fs::path outfile1a = fs::current_path() / "Lorenz_Problem_1a.csv";

    std::cout << "Initial Values & constants" << std::endl;
    std::cout << "x0 = " << x0 << " y0 = " << y0 << " z0 = " << z0 << std::endl;
    std::cout << "s = " << s << " r = " << r << " b = " << b << std::endl;
    std::cout << "t0 = " << t0 << " t final = " << tfinal << " t size = " << h << std::endl;

    LorenzRK4 lorenz_y(x0, y0, z0, t0, tfinal, s, r, b, h);
    lorenz_y.ODE_Solver();

    lorenz_y.write_csv(outfile1a);

    // std::cout << "Problem 1 printed to csv: use HW3_plots.ipynb for graphs" << std::endl << std::endl;

}


void problem2() {
    double v0, L, N, h;

    std::cout << "\nProblem 2\n" << std::endl;

    v0 = (50)*Q;         // puts V0 into eV
    N = 1000;            // number of points in the well
    h = L / N;           // step size of the well
    // other constants defined in Problem2.h file

    // fs::path testout = fs::current_path() / "TestWave.csv";
    
    // ZeroWell test(v0, 0, 5.2918*a, N, 0, Q);
    // test.secant_ODE(Q/1000);
    // double test_e = test.get_E2() / Q;
    // std::cout << " E/e = " << test_e << "ev" << std::endl;
    // test.write_function(testout);


    std::cout << "Harmonic Oscillator" << std::endl;
    HarmonicWell harmg(v0, -10*a, 10*a, N, 0.0, Q);
    harmg.secant_ODE(Q/1000.0);
    double harmg_e = harmg.get_E2() / Q;
    std::cout << "E g state = " << harmg_e << " ev" << std::endl;
    harmg.write_function(fs::current_path() / "Harmonic_Ground.csv");

    std::cout << std::endl;

    HarmonicWell harm1(v0, -10*a, 10*a, N, harmg_e*Q, 3*harmg_e*Q);
    harm1.secant_ODE(Q/1000.0);
    double harm1_e = harm1.get_E2() / Q;
    std::cout << "E 1 state = " << harm1_e << " ev" << std::endl;
    harm1.write_function(fs::current_path() / "Harmonic_First.csv");

    std::cout << std::endl;

    HarmonicWell harm2(v0, -10*a, 10*a, N, 4*harmg_e*Q, 7*harmg_e*Q);
    harm2.secant_ODE(Q/1000.0);
    double harm2_e = harm2.get_E2() / Q;
    std::cout << "E 2 state = " << harm2_e << " ev" << std::endl;
    harm2.write_function(fs::current_path() / "Harmonic_Second.csv");

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Anharmonic Oscillator" << std::endl;
    AnharmonicWell anharmg(v0, -10*a, 10*a, N, 0.0, Q);
    anharmg.secant_ODE(Q/1000.0);
    double anharmg_e = anharmg.get_E2() / Q;
    std::cout << "E g state = " << anharmg_e << " ev" << std::endl;
    anharmg.write_function(fs::current_path() / "Anharmonic_Ground.csv");

    std::cout << std::endl;

    AnharmonicWell anharm1(v0, -10*a, 10*a, N, 400*Q, 600*Q);
    anharm1.secant_ODE(Q/1000.0);
    double anharm1_e = anharm1.get_E2() / Q;
    std::cout << "E 1 state = " << anharm1_e << " ev" << std::endl;
    anharm1.write_function(fs::current_path() / "Anharmonic_First.csv");

    std::cout << std::endl;

    AnharmonicWell anhamr2(v0, -10*a, 10*a, N, 600*Q, 1000*Q);
    anhamr2.secant_ODE(Q/1000.0);
    double anhamr2_e = anhamr2.get_E2() / Q;
    std::cout << "E 2 state = " << anhamr2_e << " ev" << std::endl;
    anhamr2.write_function(fs::current_path() / "Anharmonic_Second.csv");

    std::cout << std::endl;
}


void problem3() {
    std::vector<std::vector<double>> A, A_prime;
    std::vector<double> v, voltages;

    std::cout << "Problem 3\n" << std::endl;

    std::cout << "N = 6 resistors" << std::endl;
    std::tie(A, v) = init_resistors(6, 5.0);

    A_prime = A2A(A);

    voltages = banded_matrix(A_prime, v, 2, 2);

    for (int i = 0; i < voltages.size(); i++) {
        std::cout << "V" << i + 1 << " = " << voltages[i] << std::endl;
    }

    std::cout << "Because of the unfesability of outputing 10,000 voltages to the console those will be output to a seperate .csv file" << std::endl;

    std::cout << "N = 10,000 resistors" << std::endl;
    std::tie(A, v) = init_resistors(10000, 5.0);

    A_prime = A2A(A);

    voltages = banded_matrix(A_prime, v, 2, 2);
    write_voltages(fs::current_path() / "10000_Voltages.csv", voltages);


}


int main() {
    std::cout << "Ethan Speakman HW3\n please note that the contents of the C++ file are output to .csv\nGraphs are made in the .ipynb file using Python" << std::endl << std::endl;
    // std::cout << "hbar = " << hbar << std::endl;

    // std::cout << "Problem 1 turned off right now" << std::endl;
    problem1();
    // std::cout << "Problem 2 turned off right now" << std::endl;
    problem2();
    // std::cout << "Doing problem 3" << std::endl;
    problem3();


}