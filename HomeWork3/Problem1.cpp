#include "Problem1.h"


LorenzRK4::LorenzRK4(double x0, double y0, double z0, double t0, double tfinal, double s, double r, double b, double h): 
    s(s), r(r), b(b), h(h), tfinal(tfinal) {
        // std::cout << "Init" << std::endl;
        tmax = tfinal / h;
        tmax = (int)tmax;
        // std::cout << tfinal << ", " << h << ", " << tmax << std::endl;
        // xn.resize(tmax), yn.resize(tmax), zn.resize(tmax), tn.resize(tmax);
        xn.push_back(x0), yn.push_back(y0), zn.push_back(z0); tn.push_back(t0);
        // std::cout << "Finished init" << std::endl;
        }


LorenzRK4::~LorenzRK4() {}


void LorenzRK4::ODE_Solver() {

    // std::cout << "Starting for loop with tmax = " << tmax << std::endl;
    // std::cout << "s = " << s << " r = " << r << " b = " << b << std::endl;
    for (int t = 0; t < (int)tmax; t++) {
        // std::cout << "time step " << t << std::endl;
        time_step(t);
    }
}


void LorenzRK4::time_step(int t) {
    std::array<double, 4> k_table, l_table, m_table;
    double k1, l1, m1;
    double k2, l2, m2;
    double k3, l3, m3;
    double k4, l4, m4;
    // double xnn, ynn, znn;

    // std::cout << "\t" << "xn = " << xn[t] << " yn = " << yn[t] << " zn = " << zn[t] << std::endl;

    // yes this looks odd but I was trying to find a math error and lining things up helped a lot
    k1 = dxdt(xn[t]               , yn[t]               , zn[t]);
    l1 = dydt(xn[t]               , yn[t]               , zn[t]);
    m1 = dzdt(xn[t]               , yn[t]               , zn[t]);

    // std::cout << "k1 = " << k1 << " l1 = " << l1 << " m1 = " << m1 <<std::endl;

    k2 = dxdt(xn[t] + (k1 * h / 2), yn[t] + (l1 * h / 2), zn[t] + (m1 * h / 2));
    l2 = dydt(xn[t] + (k1 * h / 2), yn[t] + (l1 * h / 2), zn[t] + (m1 * h / 2));
    m2 = dzdt(xn[t] + (k1 * h / 2), yn[t] + (l1 * h / 2), zn[t] + (m1 * h / 2));

    // std::cout << "k2 = " << k2 << " l2 = " << l2 << " m2 = " << m2 <<std::endl;

    k3 = dxdt(xn[t] + (k2 * h / 2), yn[t] + (l2 * h / 2), zn[t] + (m2 * h / 2));
    l3 = dydt(xn[t] + (k2 * h / 2), yn[t] + (l2 * h / 2), zn[t] + (m2 * h / 2));
    m3 = dzdt(xn[t] + (k2 * h / 2), yn[t] + (l2 * h / 2), zn[t] + (m2 * h / 2));

    // std::cout << "k3 = " << k3 << " l3 = " << l3 << " m3 = " << m3 <<std::endl;

    k4 = dxdt(xn[t] + (k3 * h), yn[t] + (l3 * h), zn[t] + (m3 * h));
    l4 = dydt(xn[t] + (k3 * h), yn[t] + (l3 * h), zn[t] + (m3 * h));
    m4 = dzdt(xn[t] + (k3 * h), yn[t] + (l3 * h), zn[t] + (m3 * h));

    // std::cout << "k4 = " << k4 << " l4 = " << l4 << " m4 = " << m4 <<std::endl;

    xn.push_back(xn[t] + ((h / 6) * (k1 + (2*k2) + (2*k3) + k4)));
    yn.push_back(yn[t] + ((h / 6) * (l1 + (2*l2) + (2*l3) + l4)));
    zn.push_back(zn[t] + ((h / 6) * (m1 + (2*m2) + (2*m3) + m4)));
    tn.push_back(tn[t] + h);

    // std::cout << "xn+1 = " << xn[t + 1] << " yn+1 = " << yn[t + 1] << " zn+1 = " << zn[t + 1] << std::endl;

    // std::cout << "\t" << xn[t] << ", " << yn[t] << ", " << zn[t] << std::endl;
    // std::cout << "\t\t" << xn[t + 1] << ", " << yn[t + 1] << ", " << zn[t + 1] << std::endl;

}


double LorenzRK4::dxdt(double x, double y, double z) {return (s*(y - x));}

double LorenzRK4::dydt(double x, double y, double z) {return (x*(r - z) - y);}

double LorenzRK4::dzdt(double x, double y, double z) {return ((x*y) - (b*z));}


void LorenzRK4::time(int t) {tn[t] = tn[t - 1] + h;};


void LorenzRK4::init_table() {
    blank_table.resize(3);
    blank_table[X].resize(4), blank_table[Y].resize(4), blank_table[Z].resize(4);
}


void LorenzRK4::write_csv(fs::path outfile) {

    std::cout << "Writing to file:\n" << outfile << std::endl;
    // std::cout << "tmax = " << tmax << std::endl;

    int cols = 4;

    std::ofstream fileout(outfile);

    std::string header = "t,x,y,z\n";
    fileout << header;

    for (int t = 0; t < tmax; t++) {

        std::string row_data;
        row_data = std::to_string(tn[t]) + "," + std::to_string(xn[t]) + "," + std::to_string(yn[t]) + "," + std::to_string(zn[t]) + "\n";

        fileout << row_data;
    }

    fileout.close();
    std::cout << "Finished writing to\n" << outfile << std::endl; 

}
