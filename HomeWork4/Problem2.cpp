#include "Problem2.h"


IsingLattice::IsingLattice(int N, double T, double J, double h, int seed): N(N), T(T), J(J), h(h), seed(seed) {
	// std::cout << "Initalizing Lattice for lattice\n\tN =  " << N << " T = " << T << " J = " << J << " h = " << h << std::endl << std::endl;
	init_lattice();
	// print_values();
	// print_lattice();

    // double temp = IsingLattice::T;

	// IsingLattice::T = 10;
	// // std::cout << "Warming up Lattice with T = " << IsingLattice::T << std::endl;

	// for (int t = T; t >= temp; t--){
	// 	init_holtzman();
	// 	init_joltzman();
	// 	// print_lattice();
	// }
	// IsingLattice::T = temp;
	// // std::cout << "Finished Warming up Lattice" << std::endl;
	// // std::cout << "Lattice size = " << N << "\nTemp = " << IsingLattice::T << "\nJ = " << J << " h = " << h << std::endl;

	analytic();
	// std::cout << "Analytic Solution = " << get_analytic() << std::endl;
	init_holtzman();
    init_joltzman();
	// print_energies();

	for (int s = 0; s < 10000; s++) {sweep();}
	// print_values();
	// print_lattice();

};



IsingLattice::~IsingLattice(){};



void IsingLattice::init_holtzman() {
    for (int i = -1; i <=1; i++ ) {
        holtzman[i] = exp((-2 * h * i) / T);
    }
}


void IsingLattice::init_joltzman() {
    for (int i = -2; i <= 2; i++) {
        joltzman[2*i] = exp((-2 * J * i) / T);
    }
}


void IsingLattice::sweep() {
	for (int n = 0; n < N; n++) {
		int spin_index = random_int(0, N, seed = seed);
		// std::cout << "\tspin index = " << spin_index << "/" << N << " initial spin = " << lattice[spin_index] << std::endl;
		int spin = delta_energy(spin_index);
		// std::cout << "\tspin key = " << spin << std::endl;
		double e = get_boltzman(spin, spin_index);
		// std::cout << "\tEnergy Probability = " << e << " spin index = " << spin_index << "/" << N << std::endl;
		int flip = Metropolis(e);
		// std::cout << "\tflip = " << flip << " spin index = " << spin_index << "/" << N << std::endl;
		flip_spin(flip, spin_index);
		// std::cout << "\tflip = " << flip << " spin index = " << spin_index << "/" << n << " final spin = " << lattice[spin_index] << std::endl << std::endl;
	}
}



int IsingLattice::delta_energy(int spin) {
	int e;

	int left = (spin + 1) % N;      // left one
	int right = (spin - 1 + N) % N; // right one

	e = (lattice[left] + lattice[right]) * lattice[spin]; 

	return e;
}



double IsingLattice::get_boltzman(int energy, int spin_index) {
	double delta_e;

	// std::cout << "\t\tKeys = " << 2*energy << " " << lattice[spin_index] << std::endl;

	delta_e = joltzman[2*energy] * holtzman[lattice[spin_index]];
	// std::cout << "\t\t" << joltzman[2*energy] << " " << holtzman[lattice[spin_index]] << std::endl;
	
	return delta_e;

}



void IsingLattice::flip_spin(int flip, int spin) {
	if (flip < 0) {
		// std::cout << "\t\tMove Accepted" << std::endl;
		lattice[spin] = lattice[spin] * flip;
		Mag += 2*lattice[spin];
		mag = Mag / (double)N;
	}
	
	// std::cout << "\t\t M = " << Mag << " mag = " << mag << std::endl;
}



int IsingLattice::Metropolis(float delta_e)
{
	int flip;
    float roll = ran2(seed);
	// std::cout << "\t\t r value = " << roll << std::endl;
    if (roll <= delta_e) {flip = -1;}
	else {flip = 1;}

    return flip;

}



void IsingLattice::init_magnitization() {
	Mag = 0;

	for (int n = 0; n < N; n++) {
			Mag += lattice[n];
	}

	mag = Mag / (double)N;
}



void IsingLattice::analytic() {
	double numerator = sinh(h / T);
	double denominator = sqrt(pow(sinh(h / T), 2) + exp(-4 * J / T));

	mag_analytic = numerator / denominator;
}



double IsingLattice::get_analytic() {
	return mag_analytic;
}



void IsingLattice::init_lattice() {
	lattice.resize(N);

	for (int n = 0; n < N; n++) {
		int spin = random_int(0, 2, seed);
		if (spin < 1) {spin = -1;}
		lattice[n] = spin;
	}

	init_magnitization();
}



std::tuple<double, double> IsingLattice::get_Mag() {
	return {Mag, mag};
}



void IsingLattice::print_lattice() {
	std::cout << "Printing Lattice" << std::endl;
	for (int n = 0; n < N; n++) {
		if (lattice[n] > 0) {std::cout << "+" << std::endl;}
		else {std::cout << "-" << std::endl;}
	}
}



void IsingLattice::print_values() {
	std::cout << "Lattie with:\n\tN = " << N << " T = " << T << " J = " << J << " h = " << h << " seed = " << seed << std::endl;
}


void IsingLattice::print_energies() {
	std::cout << "Holtzman Energy:" << std::endl;
	for (auto const [energy, prob]: holtzman) {std::cout << energy << ": " << prob << std::endl;}
	
	std::cout << "Joltzman Energy" << std::endl;
	for (auto const [energy, prob]: joltzman) {std::cout << energy << ": " << prob << std::endl;}

}



