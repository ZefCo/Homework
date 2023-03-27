#include "Problem2.h"

ZeroWell::ZeroWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
    V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
    potential_points();
};


ZeroWell::~ZeroWell(){};


double ZeroWell::V(double x, double a) {
    return 0.0;
}


std::tuple<double, double> ZeroWell::eigen_energy(double e) {
    double psi0, phi0;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi0 = 0.0;
    phi0 = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
        std::tie(k_psi_1, k_phi_1) = shoot(psi0, phi0, well_points[i], e);
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_1*0.5), phi0 + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_2*0.5), phi0 + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_3*0.5), phi0 + (k_phi_3*0.5), well_points[i] + (0.5*h), e);        

        psi0 += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_1)/6;
        phi0 += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_1)/6;
    }

    return {psi0, phi0};

}


void ZeroWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> ZeroWell::shoot(double psi, double phi, double x, double e) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e);

    return {fpsi, fphi};
}


void ZeroWell::secant_ODE(double target) {
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = eigen_energy(E1);

    int safty = 0;

    while (abs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = eigen_energy(E2);

        dummy1 = E2;
        E2 = E2 - psi2*(E2 - E1)/(psi2 - psi1);
        E1 = dummy1;

        safty += 1;
        if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    } 

}


double ZeroWell::get_E2() {return E2;}




HarmonicWell::HarmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
    V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
    potential_points();
};


HarmonicWell::~HarmonicWell(){};


double HarmonicWell::V(double x, double a) {
    return V0 * (x*x)/(a*a);
}


std::tuple<double, double> HarmonicWell::eigen_energy(double e) {
    double psi0, phi0;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi0 = 0.0;
    phi0 = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
        std::tie(k_psi_1, k_phi_1) = shoot(psi0, phi0, well_points[i], e);
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_1*0.5), phi0 + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_2*0.5), phi0 + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_3*0.5), phi0 + (k_phi_3*0.5), well_points[i] + (0.5*h), e);        

        psi0 += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_1)/6;
        phi0 += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_1)/6;
    }

    return {psi0, phi0};

}


void HarmonicWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> HarmonicWell::shoot(double psi, double phi, double x, double e) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e);

    return {fpsi, fphi};
}


void HarmonicWell::secant_ODE(double target) {
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = eigen_energy(E1);

    int safty = 0;

    while (abs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = eigen_energy(E2);

        dummy1 = E2;
        E2 = E2 - psi2*(E2 - E1)/(psi2 - psi1);
        E1 = dummy1;

        safty += 1;
        if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    } 

}


double HarmonicWell::get_E2() {return E2;}



AnharmonicWell::AnharmonicWell(double V0, double Ll, double Lr, double N, double E1, double E2): 
    V0(V0), Ll(Ll), Lr(Lr), N(N), E1(E1), E2(E2), h((Lr - Ll) / N) {
    potential_points();
};


AnharmonicWell::~AnharmonicWell(){};


double AnharmonicWell::V(double x, double a) {
    return V0 * (x*x*x*x)/(a*a*a*a);
}


std::tuple<double, double> AnharmonicWell::eigen_energy(double e) {
    double psi0, phi0;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi0 = 0.0;
    phi0 = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
        std::tie(k_psi_1, k_phi_1) = shoot(psi0, phi0, well_points[i], e);
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_1*0.5), phi0 + (k_phi_1*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_2*0.5), phi0 + (k_phi_2*0.5), well_points[i] + (0.5*h), e);        
        std::tie(k_psi_2, k_phi_3) = shoot(psi0 + (k_phi_3*0.5), phi0 + (k_phi_3*0.5), well_points[i] + (0.5*h), e);        

        psi0 += (k_psi_1 + (2*k_psi_2) + (2*k_psi_3) + k_psi_1)/6;
        phi0 += (k_phi_1 + (2*k_phi_2) + (2*k_phi_3) + k_phi_1)/6;
    }

    return {psi0, phi0};

}


void AnharmonicWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> AnharmonicWell::shoot(double psi, double phi, double x, double e) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e);

    return {fpsi, fphi};
}


void AnharmonicWell::secant_ODE(double target) {
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = eigen_energy(E1);

    int safty = 0;

    while (abs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = eigen_energy(E2);

        dummy1 = E2;
        E2 = E2 - psi2*(E2 - E1)/(psi2 - psi1);
        E1 = dummy1;

        safty += 1;
        if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    } 

}


double AnharmonicWell::get_E2() {return E2;}
