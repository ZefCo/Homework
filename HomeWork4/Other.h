#include <vector>
#include <math.h>
#include <ctime>
#include <chrono>
#include <iostream>


double ran2(int &idum);

int random_int(int min, int max, int seed);

double percent_diff(double a, double b);

int seed_portion(int OoM);

// Generates a seed to put into ran2. I don't care what the seed is so I just pull from the system clock
int gen_seed();

double ensemble(std::vector<double> a);