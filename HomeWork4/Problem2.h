#pragma once
#include <vector>
#include <map>
#include <tuple>
#include <math.h>
#include <iostream>
#include "Other.h"


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



// // Taken from Numerical Recipies. Adjusted so I don't have to import every single NR header file.
// double ran2(int &idum);


// // Generate random int: by default this is from [0, 1]
// // Note it actually generates a number between [0, 2) so when being used it should thought of as [min, max + 1)
// int random_int(int min, int max, int seed);