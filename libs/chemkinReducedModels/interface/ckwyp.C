#include "ckwyp.H"

void ChemkinFormationRates(const int ns, const double T, const double P, double* Y, double* R)
{
	// Workspace arrays
	int ICKWRK[500];
	double RCKWRK[500];
	
	// Conversion of pressure from Pa to dynes/cm2
	double p = P*10.;

	// Conversion of temperature
	double t = T; 

	#if defined(_WIN32) || defined(_WIN64) 

	CKWYP(&p, &t, Y, ICKWRK, RCKWRK, R);

	#else

	ckwyp_(&p, &t, Y, ICKWRK, RCKWRK, R);

	#endif

	// Conversion of molar formation rates from mol/cm3/s to kmol/m3/s
	for (int i=0;i<ns;i++)
		R[i] *= 1000.;
}
