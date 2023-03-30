#include "Problem2.h"

ZeroWell::ZeroWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
    V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
    potential_points();
};


ZeroWell::~ZeroWell(){};


double ZeroWell::V(double x, double a) {
    return 0.0;
    // return V0 * (x / Lr) * ((x / Lr) - 1);
}


std::vector<double> ZeroWell::wavefunction(double e) {
    std::vector<double> psi_func;
    double psi, phi;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi = 0.0;
    phi = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
        psi_func.push_back(psi);

        std::tie(k_psi_1, k_phi_1) = func(psi                , phi                , well_points[i]          , e);
        k_psi_1 = h * k_psi_1, k_phi_1 = h * k_phi_1;
        std::tie(k_psi_2, k_phi_2) = func(psi + (k_psi_1*0.5), phi + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
        k_psi_2 = h * k_psi_2, k_phi_2 = h * k_phi_2;
        std::tie(k_psi_3, k_phi_3) = func(psi + (k_psi_2*0.5), phi + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
        k_psi_3 = h * k_psi_3, k_phi_3 = h * k_phi_3;
        std::tie(k_psi_4, k_phi_4) = func(psi +  k_psi_3     , phi +  k_phi_3     , well_points[i] +      h , e);
        k_psi_4 = h * k_psi_4, k_phi_4 = h * k_phi_4;

        // not sure why but I have to multiply by h first... else it blows up to inf and everything fails

        psi += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_4)/6;
        phi += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_4)/6;
    }

    return psi_func;

}


void ZeroWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> ZeroWell::func(double psi, double phi, double x, double E) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - E) * psi;

    return {fpsi, fphi};
}


void ZeroWell::secant_ODE(double target) {
    std::vector<double> psi1_func, psi2_func;
    double psi1, psi2;
    double dummy1, dummy2;

    psi2_func = wavefunction(E1);
    psi2 = psi2_func.back();

    int safty = 0;

    while (fabs(E1 - E2) > target) {
        psi1 = psi2;

        psi2_func = wavefunction(E2);
        psi2 = psi2_func.back();
        // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

        dummy1 = E2;
        E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
        E1 = dummy1;

        safty += 1;
        if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    }

    unno_psi = psi2_func;
    normalize();
}


double ZeroWell::get_E2() {return E2;}


void ZeroWell::normalize() {
    std::vector<double> psipsi;
    double s;
    double dummy = unno_psi.size() / 2;  // doing it like this in case I put an odd number of well points sometime in the future
    int n = (int)dummy;

    for (int i = 0; i < n; i++) {
        psipsi.push_back(unno_psi[i] * unno_psi[i]);
    }

    s = psipsi[0] + psipsi[n];

    for (int i = 1; i < n/2; i++) {
        s += 4 * psipsi[2*i - 1];
    }

    for (int i = 1; i < (n/2) - 1; i++) {
        s += 2 * psipsi[2*i];
    }

    s = 2 * sqrt((1 / 3) * s);

    for (int i = 0; i < unno_psi.size(); i++) {
        norm_psi[i] = unno_psi[i] / s;
    }


}


void ZeroWell::write_function(fs::path outfile) {
    std::cout << "Writing to file:\n" << outfile << std::endl;
    // std::cout << "tmax = " << tmax << std::endl;

    int cols = 2;

    std::ofstream fileout(outfile);

    std::string header = "X,Phi\n";
    fileout << header;

    for (int x = 0; x < norm_psi.size(); x++) {

        std::string row_data;
        row_data = std::to_string(well_points[x]) + "," + std::to_string(norm_psi[x]) + "\n";

        fileout << row_data;
    }

    fileout.close();
    std::cout << "Finished writing to\n" << outfile << std::endl; 
}










// HarmonicWell::HarmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
//     V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
//     potential_points();
// };


// HarmonicWell::~HarmonicWell(){};


// double HarmonicWell::V(double x, double a) {
//     return V0 * (x*x) / (a*a);
//     // return V0 * (x / Lr) * ((x / Lr) - 1);
// }


// std::vector<double> HarmonicWell::wavefunction(double e) {
//     double psi, phi;
//     double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
//     double k_phi_1, k_phi_2, k_phi_3, k_phi_4;

//     std::vector<double> p_function(well_points.size());

//     psi = 0.0;
//     phi = 1.0;

//     for (int i = 0; i < well_points.size(); i++) {
//         p_function[i] = psi;
//         std::tie(k_psi_1, k_phi_1) = func(psi                , phi                , well_points[i]          , e);
//         k_psi_1 = h * k_psi_1, k_phi_1 = h * k_phi_1;
//         std::tie(k_psi_2, k_phi_2) = func(psi + (k_psi_1*0.5), phi + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
//         k_psi_2 = h * k_psi_2, k_phi_2 = h * k_phi_2;
//         std::tie(k_psi_3, k_phi_3) = func(psi + (k_psi_2*0.5), phi + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
//         k_psi_3 = h * k_psi_3, k_phi_3 = h * k_phi_3;
//         std::tie(k_psi_4, k_phi_4) = func(psi +  k_psi_3     , phi +  k_phi_3     , well_points[i] +      h , e);
//         k_psi_4 = h * k_psi_4, k_phi_4 = h * k_phi_4;

//         // not sure why but I have to multiply by h first... else it blows up to inf and everything fails

//         psi += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_4)/6;
//         phi += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_4)/6;
//     }

//     return p_function;

// }


// void HarmonicWell::potential_points() {
//     for (double i = Ll; i < Lr; i+= h) {
//         well_points.push_back(i);
//     }
// }


// std::tuple<double, double> HarmonicWell::func(double psi, double phi, double x, double e) {

//     double fpsi, fphi;

//     fpsi = phi;

//     fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e) * psi;

//     return {fpsi, fphi};
// }


// void HarmonicWell::secant_ODE(double target) {
//     std::vector<double> psi_func_1;
//     double psi1, psi2;
//     double dummy1, dummy2;

//     psi_func_1 = wavefunction(E1);
//     psi1 = psi_func_1.back();

//     int safty = 0;

//     while (fabs(E1 - E2) > target) {
//         psi1 = psi2;

//         std::tie(psi2, dummy2) = wavefunction(E2);
//         // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

//         dummy1 = E2;
//         E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
//         E1 = dummy1;

//         safty += 1;
//         // if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
//     } 
//     std::cout << "Iterations = " << safty << std::endl;
// }


// double HarmonicWell::get_E2() {return E2;}










// AnharmonicWell::AnharmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
//     V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
//     potential_points();
// };


// AnharmonicWell::~AnharmonicWell(){};


// double AnharmonicWell::V(double x, double a) {
//     return V0 * (x*x*x*x)/(a*a*a*a);
// }


// std::tuple<double, double> AnharmonicWell::wavefunction(double e) {
//     double psi, phi;
//     double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
//     double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


//     psi = 0.0;
//     phi = 1.0;

//     for (int i = 0; i < well_points.size(); i++) {
//         std::tie(k_psi_1, k_phi_1) = func(psi                , phi                , well_points[i]          , e);
//         k_psi_1 = h * k_psi_1, k_phi_1 = h * k_phi_1;
//         std::tie(k_psi_2, k_phi_2) = func(psi + (k_psi_1*0.5), phi + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
//         k_psi_2 = h * k_psi_2, k_phi_2 = h * k_phi_2;
//         std::tie(k_psi_3, k_phi_3) = func(psi + (k_psi_2*0.5), phi + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
//         k_psi_3 = h * k_psi_3, k_phi_3 = h * k_phi_3;
//         std::tie(k_psi_4, k_phi_4) = func(psi +  k_psi_3     , phi +  k_phi_3     , well_points[i] +      h , e);
//         k_psi_4 = h * k_psi_4, k_phi_4 = h * k_phi_4;

//         // not sure why but I have to multiply by h first... else it blows up to inf and everything fails

//         psi += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_4)/6;
//         phi += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_4)/6;
//     }

//     return {psi, phi};

// }


// void AnharmonicWell::potential_points() {
//     for (double i = Ll; i < Lr; i+= h) {
//         well_points.push_back(i);
//     }
// }


// std::tuple<double, double> AnharmonicWell::func(double psi, double phi, double x, double e) {

//     double fpsi, fphi;

//     fpsi = phi;

//     fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e) * psi;

//     return {fpsi, fphi};
// }


// void AnharmonicWell::secant_ODE(double target) {
//     double psi1, psi2;
//     double dummy1, dummy2;

//     std::tie(psi2, dummy2) = wavefunction(E1);

//     int safty = 0;

//     while (fabs(E1 - E2) > target) {
//         psi1 = psi2;

//         std::tie(psi2, dummy2) = wavefunction(E2);
//         // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

//         dummy1 = E2;
//         E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
//         E1 = dummy1;

//         safty += 1;
//         if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
//     } 
// }

// double AnharmonicWell::get_E2() {return E2;}
