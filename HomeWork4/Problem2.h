#include <vector>
#include <map>

int seed;


class IsingLattice {
    public:
        IsingLattice(double N, double T, double J, double h);
        ~IsingLattice();

        std::map<int, double> joltzman;
        std::map<int, double> holtzman;


    private:
        // Coupeling Constant
        double J;
        // Tempurature
        double T;
        // external feild
        double h;
        // size of lattice
        double N;

        void init_joltzman();

        // Initialize the arrays that hold the Boltzman values
        // holtzman for H Boltzman (not a Dune reference... or is it?)
        void init_holtzman();

};



// Taken from Numerical Recipies. Adjusted so I don't have to import every single NR header file.
double ran2(int &idum);


// Generate random int: by default this is from [0, 1]
// Note it actually generates a number between [0, 2) so when being used it should thought of as [min, max + 1)
int random_int(int min = 0, int max = 2);