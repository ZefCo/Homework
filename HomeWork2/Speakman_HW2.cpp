# include <iostream>
# include <vector>
# include <random>
# include <math.h>
# include <map>
# include <fstream>
# include <numeric>
# include <string>
# include <functional>
// # include "matplotlibcpp.h"



double bisectional_method(double x, double a, double b, double tolerence, double (*func)(double w, double x)) {

    double c;
    double fa, fb, fc;

    fa = func(a, x);
    fb = func(b, x);

    if ((fa * fb) > 0.0) {return false;}

    while ((b - a) > tolerence) {
        c = (a + b) / 2.0;

        fc = func(c, x);

        if ((fc * fa) < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return ((a + b) / 2.0);

}



int problem1() {
    // W*exp(W) = x
    // W*exp(W) - x = 0
    double tolerence = 10^-10;

    auto lambertW = [](double w, double x) -> double {return w*exp(w) - x;};




}


float jacobian() {
    return 0.0;
}


// Lecture 5 - Newton Method of Root Finding
// Then put in values for V1 and V2 and find I for Diode, then use that to find
// Voltage. Need Resistance/Inductance of Diode.
void problem2(float x0, float accuracy, int max_iterations) {
    // 

    float x = x0;

    for (int i = 0; i < max_iterations; i++) {

    };

}


int main() {

}