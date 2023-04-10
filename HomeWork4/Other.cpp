#include "Other.h"

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


int random_int(int min = 0, int max = 2, int seed = 20170520)
{

    double rval = ran2(seed);
    int i = min + ((max - min) * rval);

    return i;
}


double percent_diff(double a, double b) {
    double c;

    c = abs((a - b) / (0.5 * (a + b)));

    return c*100;
}





int seed_portion(int OoM = 1) {
    // std::chrono::duration<double> timeX, timeY;

    auto timeX = std::chrono::system_clock::now();  // autotype because whatever the fuck this is, typeid().name() seems to hate it.
    auto timeY = std::chrono::system_clock::now();

    int deltaT = (timeY - timeX).count() * (int)OoM;

    return deltaT;


}


int gen_seed() {

    int seedy1 = seed_portion();
    int seedy2 = seed_portion(1000);
    int seedy3 = seed_portion(1000000);

    return (-1)*(seedy1 + seedy2 + seedy3);
}



double ensemble(std::vector<double> a) {

    double x = 0;

    double X = a.size();

    for (int i = 0; i < a.size(); i++) {
        x += a[i];
    }

    return x / X;

}
