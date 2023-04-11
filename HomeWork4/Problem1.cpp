#include "Problem1.h"


EarthCrust::EarthCrust(double years, double TEarthCore, double TEarthSurface, int N, int L, double h): 
TEarthCore(TEarthCore), TEarthSurface(TEarthSurface), N(N), L(L), h(h), years(years) {
    init_grid();

    a = (double)L / (double)N;
    c = h * D / (a*a);
    epsilon = h / 1000;

    warmup = (years - 1.0)*tau;

    first_quarter = (years - 1.0 + 0.25) * tau;
    second_quarter = (years - 1.0 + 0.5) * tau;
    third_quarter = (years - 1.0 + 0.75) * tau;
    fourth_quarter = years*tau;
    
    tend = fourth_quarter + epsilon;

};


EarthCrust::~EarthCrust() {};


double EarthCrust::T(double t) {
    return A + B*sin(2 * PI * t / tau);
};


void EarthCrust::init_grid() {
    depth.resize(N + 1);
    depth_prime.resize(N + 1);
    
    for (int n = 0; n < N + 1; n++) {
        depth[n] = TEarthSurface;
        depth_prime[n] = TEarthSurface;
    }

    depth[0] = TEarthCore;
    depth_prime[0] = TEarthCore;
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> EarthCrust::run_sim() {
    double t = 0;    
    std::vector<double> fiq, seq, thq, foq;

    while (t < tend) {

        depth[0] = T(t);
        // std::cout << depth[0] << std::endl;
        for (int i = 1; i < N; i++) {
            depth_prime[i] = depth[i] + c*(depth[i + 1] + depth[i - 1] - 2.0*depth[i]);
        }

        std::vector<double> dummy = depth_prime;
        depth_prime = depth;
        depth = dummy;

        if (abs(t - first_quarter) < epsilon) {fiq = depth;}
        if (abs(t - second_quarter) < epsilon) {seq = depth;}
        if (abs(t - third_quarter) < epsilon) {thq = depth;}
        if (abs(t - fourth_quarter) < epsilon) {foq = depth;}

        t += h;
    }

    return {fiq, seq, thq, foq};
}


void write_problem1(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d) {
    fs::path outfile = fs::current_path() / "Problem1Data.csv";
    std::cout << "Writing to file:\n" << outfile << std::endl;

    std::ofstream fileout(outfile);

    fileout << "n";
    for (int q = 1; q < 5; q++) {fileout << ",Q_" << q;}
    fileout << "\n";

    for (int i = 0; i < a.size(); i++) {
        fileout << i << "," << a[i] << "," << b[i] << "," << c[i] << "," << d[i] << "\n";
    }

    fileout.close();

    std::cout << "Finished writing to\n" << outfile << std::endl;

}
