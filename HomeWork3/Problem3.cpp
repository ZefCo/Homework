#include "Problem3.h"


std::tuple<std::vector<std::vector<double>>, std::vector<double>> init_resistors(int N, double V) {

    std::vector<std::vector<double>> resistors;
    std::vector<double> voltages(N);

    resistors.resize(N);
    
    for (int n = 0; n < N; n++) {resistors[n].resize(N);}

    for (int n = 1; n < N - 1; n++) {resistors[n][n] = 4;}
    
    for (int n = 0; n < N - 1; n++) {
        resistors[n][n + 1] = -1; 
        resistors[n + 1][n] = -1;}
    
    for (int n = 0; n < N - 2; n++) {
        resistors[n][n + 2] = -1;
        resistors[n + 2][n] = -1;}

    resistors[0][0] = 3, resistors[N - 1][N - 1] = 3;

    voltages[0] = V, voltages[0] - V;

    // std::cout << std::endl << "Resistor Matrix" << std::endl;
    // for (int r = 0; r < N; r++){
    //     for (int c = 0; c < N; c++) {
    //         if (resistors[r][c] >= 0) {std::cout << "+" << resistors[r][c] << " ";}
    //         else {std::cout << resistors[r][c] << " ";}
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "Voltages" << std::endl;
    // for (int v = 0; v < N; v++) {
    //     std::cout << voltages[v] << " ";
    // }
    // std::cout << std::endl << std::endl;

    return {resistors, voltages};

}


std::vector<double> banded_matrix(std::vector<std::vector<double>> Ain, std::vector<double> vin, int up, int down) {
    int N;
    double div;
    std::vector<double> va;
    std::vector<std::vector<double>> Aa;

    N = vin.size();

    va.resize(N);
    Aa.resize(N);
    for (int r = 0; r < N; r++) {Aa[r].resize(N);}

    for (int v = 0; v < vin.size(); v++) {va[v] = vin[v];}
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            Aa[r][c] = Ain[r][c];
        }
    }

    for (int m = 0; m < N; m++) {
        div = Aa[up][m];

        va[m] /= div;
        for (int k = 1; k <= down + 1; k++) {
            if (m + k < N) {va[m + k] -= Aa[up + k][m] * va[m];}
        }

        for (int i = 0; i < up; i++) {
            int j = m + up - i;
            if (j < N) {
                Aa[i][j] /= div;
                for (int k = 0; k <= down + 1; k++) {
                    Aa[i + k][j] -= Aa[up + k][m] * Aa[i][j];}
            }
        }
    }

    for (int m = N - 2; m > -1; m--) {
        for (int i = 0; i < up; i++) {
            int j = m + up - i;
            if (j < N) {va[m] -= Aa[i][j] * va[j];}
        }
    }

    return va;

}




std::vector<std::vector<double>> A2A(std::vector<std::vector<double>> inA) {
    int N = inA[0].size(); // assumes a square matrix
    int row = 0;

    std::vector<std::vector<double>> outA;
    outA.resize(N);
    for (int r = 0; r < N; r++) {outA[r].resize(N);}
    
    // iterate across first row
    for (int c = N; c > -1; c--) {
        if (abs(inA[0][c]) > 0) {       // note the abs(). Originally I used x != 0, but somehow, when checking this on other computers, it would find none zero
            std::vector<double> dummy;  // values where it should be zero, meaning it created a false diagnol row. Switching to abs(x) > 0 seems to have cleared
                                        // things up. I could have used and even wider tolerance because all the original coefficents are >= 1.
            for (int i = 0; i < N - c; i++) {
                dummy.push_back(inA[i + c][i]);
            }

            if ((N - dummy.size()) > 0) {
                for (int i = 0; i <= N - dummy.size(); i++) {
                    dummy.insert(dummy.begin(), 0);
                }
            }

            for (int i = 0; i < N; i ++) {outA[row][i] = dummy[i];}
            row += 1;
        }
    }

    // iterate across first column
    // starting at 1 because the above ends at 0
    for (int r = 1; r < N; r++) {
        if (abs(inA[r][0]) > 0) {
            std::vector<double> dummy;

            for (int i = 0; i < N - r; i++) {
                dummy.push_back(inA[i][i + r]);
            }

            if ((N - dummy.size()) > 0 ) { 
                for (int i = 0; i <= N - dummy.size(); i++) {
                    dummy.push_back(0);
                }
            }

            for (int i = 0; i < N; i++) {outA[row][i] = dummy[i];}
            row += 1;

        }
    }

    // std::cout << std::endl << "Diagnol Matrix" << std::endl;
    // for (int r = 0; r < N; r++) {
    //     for (int c = 0; c < N; c++) {
    //         if (outA[r][c] >= 0) {std::cout << "+" << outA[r][c] << " ";}
    //         else {std::cout << outA[r][c] << " ";}
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    return outA;

}


void write_voltages(fs::path outfile, std::vector<double> volts) {
    std::cout << "Writing to file:\n" << outfile << std::endl;
    // std::cout << "tmax = " << tmax << std::endl;

    int cols = 2;

    std::ofstream fileout(outfile);

    std::string header = "i,Volts\n";
    fileout << header;

    for (int v = 0; v < volts.size(); v++) {
        fileout << "V" << std::to_string(v);
        fileout << ",";
        fileout << volts[v];
        fileout << "\n";

    }

    fileout.close();
    std::cout << "Finished writing to\n" << outfile << std::endl; 

}