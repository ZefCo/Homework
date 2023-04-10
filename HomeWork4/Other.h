#include <vector>
#include <math.h>
#include <ctime>
#include <chrono>
#include <iostream>


// Taken from Numerical recipies. Only modifications are to use C++ declarations as opposed to NR declarations and have like 6 more files
double ran2(int &idum);


// Uses NR to get a random int. Default return is [0, 2) and default seed is -20170520
int random_int(int min, int max, int seed);

// Finds the % diff between two values
double percent_diff(double a, double b);

// Does the actual finding of the seed from the system clock. Called in gen_seed.
int seed_portion(int OoM);

// Generates a seed to put into ran2. I don't care what the seed is so I just pull from the system clock
int gen_seed();

// Finds the average of a vector.
double ensemble(std::vector<double> a);