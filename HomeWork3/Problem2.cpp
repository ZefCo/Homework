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


std::tuple<double, double> ZeroWell::wavefunction(double e) {
    double psi, phi;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi = 0.0;
    phi = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
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

    return {psi, phi};

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
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = wavefunction(E1);

    int safty = 0;

    while (fabs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = wavefunction(E2);
        // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

        dummy1 = E2;
        E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
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
    return V0 * (x*x) / (a*a);
    // return V0 * (x / Lr) * ((x / Lr) - 1);
}


std::tuple<double, double> HarmonicWell::wavefunction(double e) {
    double psi, phi;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi = 0.0;
    phi = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
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

    return {psi, phi};

}


void HarmonicWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> HarmonicWell::func(double psi, double phi, double x, double e) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e) * psi;

    return {fpsi, fphi};
}


void HarmonicWell::secant_ODE(double target) {
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = wavefunction(E1);

    int safty = 0;

    while (fabs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = wavefunction(E2);
        // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

        dummy1 = E2;
        E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
        E1 = dummy1;

        safty += 1;
        // if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    } 
    std::cout << "Iterations = " << safty << std::endl;
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


std::tuple<double, double> AnharmonicWell::wavefunction(double e) {
    double psi, phi;
    double k_psi_1, k_psi_2, k_psi_3, k_psi_4;
    double k_phi_1, k_phi_2, k_phi_3, k_phi_4;


    psi = 0.0;
    phi = 1.0;

    for (int i = 0; i < well_points.size(); i++) {
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

    return {psi, phi};

}


void AnharmonicWell::potential_points() {
    for (double i = Ll; i < Lr; i+= h) {
        well_points.push_back(i);
    }
}


std::tuple<double, double> AnharmonicWell::func(double psi, double phi, double x, double e) {

    double fpsi, fphi;

    fpsi = phi;

    fphi = (2 * me / (hbar*hbar)) * (V(x, a) - e) * psi;

    return {fpsi, fphi};
}


void AnharmonicWell::secant_ODE(double target) {
    double psi1, psi2;
    double dummy1, dummy2;

    std::tie(psi2, dummy2) = wavefunction(E1);

    int safty = 0;

    while (fabs(E1 - E2) > target) {
        psi1 = psi2;

        std::tie(psi2, dummy2) = wavefunction(E2);
        // std::cout << "Psi1 = " << psi1 << " Psi2 = " << psi2 << std::endl;

        dummy1 = E2;
        E2 = E2 - (psi2 * (E2 - E1) / (psi2 - psi1));
        E1 = dummy1;

        safty += 1;
        if (safty > 100) {std::cout << "safty rail hit" << std::endl; std::cout << "E1 = " << E1 << " E2 = " << E2 << std::endl; break;}
    } 
}

double AnharmonicWell::get_E2() {return E2;}
