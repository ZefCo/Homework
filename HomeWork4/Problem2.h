#pragma once
#include <vector>
#include <map>
#include <tuple>
#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "Other.h"

namespace fs = std::filesystem;



class IsingLattice {
    public:
        IsingLattice(int N, double T, double J, double h, int seed);
        ~IsingLattice();

        // Stores the joltzman values. Instead of finding these every loop I store them as a map and index as needed
        std::map<int, double> joltzman;
        // Stores the holtzman values. Instead of finding these every loop I store them as a map and index as needed
        std::map<int, double> holtzman;

        // Preforms 1 sweep across the lattice
        void sweep();

        // Finds the analytic solution
        double get_analytic();

        // gets the energy from the J and h sums. Indexes both the joltzman and holtzman values, then multiplies the values.
        double get_boltzman(int energy, int spin);

        // finds the delta of the spins. This actually returns an int, which is the key to the get_boltzman function. 
        // So with the value of the this you can index to find what the energy is.
        int delta_energy(int spin);

        // Determines if the spin will flip based on Metropolis Rates, returns a -1 or 1
        int Metropolis(float delta_e);

        // flips the spin if the flip input is -1. Yes I could probably put this in with the Metropolis function but I didn't.
        void flip_spin(int flip, int spin);

        // returns the magnitization of the lattice, both Mag and Mag/spin, as a tuple
        std::tuple<double, double> get_Mag();


        // Used these when I was debugging
        void print_lattice();
        void print_values();
        void print_energies();

    private:
        // Coupeling Constant
        double J;
        // Tempurature
        double T;
        // external feild
        double h;
        // size of lattice
        int N;

        // Initialize the arrays that hold teh Boltzman values joltzman for J Boltzman (sadly this is not a reference to anything)
        // Runs from -4 to +4, though I only need the values for -4, 0, and +4
        // It's just easier to get them all (and safer, in case I did the math wrong)
        void init_joltzman();

        // Initialize the arrays that hold the Boltzman values holtzman for H Boltzman (not a Dune reference... or is it?)
        // Runs from -2 to +2 even though I only need the values for +2 and -2
        void init_holtzman();

        // Creates and populates the lattice with spins
        void init_lattice();

        // Finds the inital Magnitization of the lattice
        void init_magnitization();

        // Finds the analytic solution for the 1D lattice
        void analytic();


        // The declaration of the lattice, initialized above.
        std::vector<int> lattice;

        // Energy for the system
        double energy;
        // magnitization per spin
        double mag;
        // Overall Magnitization
        double Mag;

        // seed for the random number generator. 
        // I was conceptually having trouble finding a place for it so I just stuck it in as an attribute. 
        // Probably not the best place for it but it works.
        int seed;

        // The analytic solution for the 1D Ising Model
        double mag_analytic;

};



// Iterates from N = 10 to N = 1000 to find a decent value for N. Will generate a number of seeds and average across those seeds to find <|m|>.
// Once it finds one within a certain trehshold it stops the loop (so you don't have to sit threw 990 iterations)
void finding_N(double T, double J, double h, double thresh, int runs);


// Does the same thing as finding_N but only for one value of N. Good for vaildating the value you found. 
// Because of the nature of the PRNG you may need to run it once or twice: using N = 75 I had a percent diff from the analytic solution of 9, 5, and 1%
// with 5 runs <2% with 10 runs. More runs is better, but longer.
void validate_N(int N, double T, double J, double h, int runs);



// Forget finding N, I was doing in a convoluted way. Just doing a bunch of sweeps for several N and ploting.
// Has one issue I'm not going to bother fixing: you can iterate outside of Max N which will cause it to crash. Think about it:
// if you choose a increment value that goes beyond Max N (like min = 10, max = 50, increment = 25), it crashes. 
// But since those values are hard coded to the main I'm not going to worry about. Maybe sometime in the future.
void finding_N_better(double T, double J, double h, int min_N, int max_N, int increment);

// Writes the Data out to a csv file. Used in both the finding_N_Better and h_values script
void write_problem2(fs::path outfile, std::vector<std::vector<double>> m, std::vector<std::string> header);


// Varries the h parameter, but takes in a T, J, N. Has the same issue that finding_N_better does
void h_values(int N, double T, double J, double min_h, double max_h, double increment);