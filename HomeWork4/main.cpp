#include <iostream>
#include "Problem2.h"
// #include "Other.h"


// g++ main.cpp Problem2.cpp Other.cpp -o SpeakmanHW4.exe
// g++ main.cpp Problem2.cpp Other.cpp -o SpeakmanHW4.out


void Problem1() {
    
}



void Problem2() {
    double T = 1;
    double J = 1;
    double h = 0.1;
    double thresh = 0.1;
    int seperate_seeds = 10;

    // finding_N(T, J, h, thresh, seperate_seeds);  
    // Value was found to be between 20-40, so taking 30 and doubling it I'd say 60 is a good number to pick for N -> inifinity

    validate_N(100, T, J, h, 10);


}



int main() {
    Problem2();

}