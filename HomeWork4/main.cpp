#include <iostream>
#include "Problem2.h"
// #include "Other.h"


// g++ main.cpp Problem2.cpp Other.cpp -o SpeakmanHW4.exe

void Problem1() {
    
}



void Problem2a() {
    double T = 1;
    double J = 1;
    double h = 0.1;
    double thresh = 1;


    // IsingLattice oneD1(10, T, J, h, -20170520);
    // std::vector<double> m_vec, M_vec;
    // double ave_m, M, m;

    // for (int i = 0; i < 4; i++){
    //     for (int s = 0; s < 100000; s++) {oneD1.sweep();}
    //     oneD1.print_lattice();
    //     std::tie(M, m) = oneD1.get_Mag();
    //     // std::cout << "M = " << M << " m = " << m << std::endl;
    //     m_vec.push_back(abs(m)); M_vec.push_back(abs(M));
    //     // std::cout << std::endl;
    //     // oneD.print_lattice();
    //     // std::cout << std::endl;
    // }

    // for (int i = 0; i < m_vec.size(); i++) {std::cout << "m" << i << " = " << m_vec[i] << " M" << i << " = " << M_vec[i] << std::endl;}
    // ave_m = ensemble(m_vec);

    // std::cout << "<m> = " << ave_m << std::endl;
    // // oneD1.print_values();
    // // oneD1.print_lattice();
    // // oneD1.print_energies();
    

    // IsingLattice oneD2(12, T, J, h, -20170520);
    // std::vector<double> m_vec2, M_vec2;
    // double ave_m2, M2, m2;

    // for (int i = 0; i < 4; i++){
    //     for (int s = 0; s < 100000; s++) {oneD2.sweep();}
    //     std::tie(M, m) = oneD2.get_Mag();
    //     // std::cout << "M = " << M << " m = " << m << std::endl;
    //     m_vec2.push_back(abs(m)); M_vec2.push_back(abs(M));
    //     // std::cout << std::endl;
    //     // oneD.print_lattice();
    //     // std::cout << std::endl;
    // }

    // for (int i = 0; i < m_vec2.size(); i++) {std::cout << "m" << i << " = " << m_vec2[i] << " M" << i << " = " << M_vec2[i] << std::endl;}
    // ave_m2 = ensemble(m_vec2);

    // std::cout << "< |m| > = " << ave_m2 << std::endl;



    for (int n = 10; n < 100; n++) {
        std::vector<double> m_ave_vec;
        double can_m, pd, analytic;
        int seperate_seeds = 5;
        std::cout << "#### N = " << n << " with " << seperate_seeds << " seperate seeds" << std::endl;
        for (int e = 0; e < seperate_seeds; e++) {
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



int main() {
    Problem2a();

}