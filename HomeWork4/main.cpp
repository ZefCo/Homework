#include <iostream>
#include <tuple>
#include "Other.h"
#include "Problem1.h"
#include "Problem2.h"


// g++ main.cpp Problem2.cpp Other.cpp -o SpeakmanHW4.exe
// g++ main.cpp Problem2.cpp Other.cpp -o SpeakmanHW4.out


void Problem1() {
    double TEarthCore = 11;
    double TEarthSurface = 10;
    double h = 1.0e-2;
    double L = 20;
    double N = 100;

    double years = 10;

    std::cout << "Problem 1\nFinding the Temp Profile of the Earth's crust" << std::endl;
    std::cout << "Using Values\n\tEarth's Surface = " << TEarthSurface << " Depth = " << L << " Temp @ 20 m = " << TEarthCore << " time step = " << h << " Grid spacing = " << N << std::endl;
    std::cout << "Running for 9 years then output each quarter of the last year" << std::endl;

    std::vector<double> first_qart, second_quart, third_quart, fourth_quart;

    EarthCrust ec(years, TEarthCore, TEarthSurface, N, L, h);
    std::tie(first_qart, second_quart, third_quart, fourth_quart) = ec.run_sim();

    write_problem1(first_qart, second_quart, third_quart, fourth_quart);
    
}



void Problem2() {
    double T = 1;
    double J = 1;
    double h0 = 0.1;

    int min_N = 500;
    int max_N = 2000;
    int N_step = 50;

    int N = 2000;
    double h_min = -0.5;
    double h_max = 0.5;
    double h_step = 0.1;
    int seperate_seeds = 10;
    
    std::cout << std::endl;
    std::cout << "Problem 2\n1D Ising Model Simulation\n" << std::endl;
    std::cout << "T = " << T << " J = " << J << " h = " << h0 << std::endl;
    std::cout << "A table will be output for different values of N\nPrevious runs of this data showed N ~ 1000 gave good results for spontaneous magnitization, so this will run from N = " << min_N << " to " << max_N << " in steps of " << N_step << std::endl;
    finding_N_better(T, J, h0, min_N, max_N, N_step);
    std::cout << std::endl << "\nUsing N = " << N << ", now a table with different values of h varrying from h = " << h_min << " to " << h_max << " in steps of " << h_step << std::endl;
    h_values(N, T, J, h_min, h_max, h_step);

}



int main() {
    Problem1();
    Problem2();

}