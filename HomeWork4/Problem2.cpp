#include "Problem2.h"


IsingLattice::IsingLattice(int N, double T, double J, double h, int seed): N(N), T(T), J(J), h(h), seed(seed) {
	init_lattice();
	analytic();
	init_holtzman();
    init_joltzman();

	// for (int s = 0; s < 10000; s++) {sweep();}

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
		int spin = delta_energy(spin_index);
		double e = get_boltzman(spin, spin_index);
		int flip = Metropolis(e);
		flip_spin(flip, spin_index);
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

	delta_e = joltzman[2*energy] * holtzman[lattice[spin_index]];
	
	return delta_e;

}



void IsingLattice::flip_spin(int flip, int spin) {
	if (flip < 0) {
		lattice[spin] = lattice[spin] * flip;
		Mag += 2*lattice[spin];
		mag = Mag / (double)N;
	}
	
}



int IsingLattice::Metropolis(float delta_e)
{
	int flip;
    float roll = ran2(seed);
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



void finding_N(double T, double J, double h, double thresh, int runs) {

	// seraching for N
    for (int n = 10; n < 100; n++) {
        std::vector<double> m_ave_vec;
        double can_m, pd, analytic;
        
		std::cout << "#### N = " << n << " with " << runs << " seperate seeds" << std::endl;
        for (int r = 0; r < runs; r++) {
            std::vector<double> m_vec, M_vec;
            int seed = gen_seed();
            // int seed = -20170520;
            double ave_m;
            double M, m;
            int measures = 10;  // OK yes I should do an autocorrelation and other things but... that's a lot of work

            IsingLattice oneD(n, T, J, h, seed);

            for (int i = 0; i < measures; i++){
                double M, m;
                for (int s = 0; s < 100000; s++) {oneD.sweep();}
                std::tie(M, m) = oneD.get_Mag();
                // std::cout << "M = " << M << " m = " << m << std::endl;
                m_vec.push_back(m); M_vec.push_back(M); 
                // std::cout << "\t|m| = " << abs(m) << std::endl;
            }

            ave_m = ensemble(m_vec);
            m_ave_vec.push_back(ave_m);
            analytic = oneD.get_analytic();
        }

        can_m = ensemble(m_ave_vec);

        pd = percent_diff(can_m, analytic);

        std::cout << "< |m| > = " << can_m << " target = " << analytic << " diff = " << pd << " N = " << n << std::endl;
        // std::cout << "Length of m_vec = " << m_vec.size() << std::endl;
        if (pd <= thresh) {break;}

    }

}



void validate_N(int N, double T, double J, double h, int runs) {
    std::vector<double> m_ave_vec;
    double can_m, pd, analytic;

    std::cout << "#### N = " << N << " with " << runs << " seperate seeds" << std::endl;

    for (int r = 0; r < runs; r++) {
        std::vector<double> m_vec, M_vec;
        int seed = gen_seed();
        double ave_m;
        double M, m;
        int measures = 10;  // OK yes I should do an autocorrelation and other things but... that's a lot of work

        IsingLattice oneD(N, T, J, h, seed);

        for (int i = 0; i < measures; i++){
            double M, m;
            for (int s = 0; s < 100000; s++) {oneD.sweep();}
            std::tie(M, m) = oneD.get_Mag();
            // std::cout << "M = " << M << " m = " << m << std::endl;
            m_vec.push_back(m); M_vec.push_back(M); 
            // std::cout << "\t|m| = " << abs(m) << std::endl;
        }

        ave_m = ensemble(m_vec);
        m_ave_vec.push_back(ave_m);
        analytic = oneD.get_analytic();
    }

    can_m = ensemble(m_ave_vec);

    pd = percent_diff(can_m, analytic);

    std::cout << "< |m| > = " << can_m << " target = " << analytic << " diff = " << pd << " N = " << N << std::endl;

}


// comes out to ~ 390 = N
void finding_N_better(double T, double J, double h, int min_N, int max_N, int increment) {
    std::vector<std::vector<double>> m_data;
    int i = 0;
    std::vector<std::string> headers;

    m_data.resize(((max_N - min_N) / increment) + 1);
    headers.resize(((max_N - min_N) / increment) + 1);
    
    for (int n = min_N; n < max_N + increment; n += increment) {
		std::cout << "#### Run " << i + 1 << " N = " << n << std::endl;
        headers[i] = "N_" + std::to_string(n);
        
        std::vector<double> m_vec;

        int seed = gen_seed();
        IsingLattice oneD(n, T, J, h, seed);
        double M, m;

        for (int i = 0; i < 10000; i++){
            oneD.sweep();
            std::tie(M, m) = oneD.get_Mag();
            m_vec.push_back(m); 
        }

        m_data[i] = m_vec; i += 1;
    }

    write_problem2(fs::current_path() / "MagData_N.csv", m_data, headers);


}



void write_problem2(fs::path outfile, std::vector<std::vector<double>> m, std::vector<std::string> headers) {

    std::cout << "Writing to file:\n" << outfile << std::endl;
    // std::cout << "tmax = " << tmax << std::endl;
    std::ofstream fileout(outfile);
    int cols = headers.size();

    fileout << "t";
    for (int c = 0; c < cols; c++) {fileout << "," << headers[c];}
    fileout << "\n";

    for (int t = 0; t < m[0].size(); t++) {
        fileout << t;
    
        for (int c = 0; c < cols; c++) {
            fileout << "," << m[c][t];
        }
    
        fileout << "\n";
    }

    fileout.close();
    std::cout << "Finished writing to\n" << outfile << std::endl; 

}



void h_values(int N, double T, double J, double min_h, double max_h, double increment) {
    std::vector<std::vector<double>> m_data;
    int i = 0;
    std::vector<std::string> headers;
    int vsize = (int)((max_h - min_h) / increment);
    m_data.resize(vsize);
    headers.resize(vsize);

    for (double h = min_h; h < max_h; h += increment) {
        if (abs(h) < 0.01) {h = 0;}
    
		std::cout << "#### Run " << i + 1 << " h = " << h << std::endl;
        headers[i] = "h_" + std::to_string(h);

        std::vector<double> m_vec;

        int seed = gen_seed();
        IsingLattice oneD(N, T, J, h, seed);
        double M, m;

        for (int i = 0; i < 10000; i++){
            oneD.sweep();
            std::tie(M, m) = oneD.get_Mag();
            m_vec.push_back(m); 
        }
        m_data[i] = m_vec; i += 1;
    }

    write_problem2(fs::current_path() / "MagData_h.csv", m_data, headers);

}