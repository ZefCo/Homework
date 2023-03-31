#include <vector>
#include <array>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <tuple>
#include <filesystem>

namespace fs = std::filesystem;


// creates a matrix of resistors and voltages.
std::tuple<std::vector<std::vector<double>>, std::vector<double> > init_resistors(int N, double V);

// A C++ adaptation of the banded.py file from pg 525 of Mark Newman's COmputational Physics.
std::vector<double> banded_matrix(std::vector<std::vector<double>> Aa, std::vector<double> va, int up, int down);

// converts the diagnols of the matrix into rows of another matrix
// this will output a square matrix, but the contents of this matrix are for the banded_matrix function
std::vector<std::vector<double>> A2A(std::vector<std::vector<double>> inA);

// writes the voltages to a csv file: mainly for the 10,000 volt problem as printing those to the console is not a smart idea
void write_voltages(fs::path outfile, std::vector<double> volts);