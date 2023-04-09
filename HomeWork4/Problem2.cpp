#include "Problem2.h"


IsingLattice::IsingLattice(double N, double T, double J, double h): N(N), T(T), J(J), h(h) {
    init_holtzman();
    init_joltzman();

};



IsingLattice::~IsingLattice(){};



void IsingLattice::init_holtzman() {
    for (int i = -1; i <=1; i++ ) {
        holtzman[i] = exp((-h * i) / T);
    }
}


void IsingLattice::init_joltzman() {
    for (int i = -2; i <= 2; i++) {
        joltzman[i] = exp((-1.0 * i) / T);
    }
}


double ran2(int &idum)
{
	const int IM1=2147483563,IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	static int idum2=123456789,iy=0;
	static std::vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (idum <= 0) {
		idum=(idum==0 ? 1 : -idum);
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


int random_int(int min = 0, int max = 2)
{

    double rval = ran2(seed);
    int i = min + ((max - min) * rval);

    return i;
}