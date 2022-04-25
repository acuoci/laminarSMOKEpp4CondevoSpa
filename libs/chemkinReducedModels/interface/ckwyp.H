// Returns the molar production rates of the species given pressure, temperature and mass fractions
// P: pressure in cgs units (dynes/cm2)
// T: temperature in K
// Y: mass fractions of species
// ICKWRK: integer workspace array, dimension at least LENICK
// LENRCK: real workspace array, dimension at least LENRCK
// WDOT: formation rates of species in mol/cm3/s

#if defined(_WIN32) || defined(_WIN64) 

	extern "C" {	void CKWYP(double* P, double* T, double* Y, int* ICKWRK, double* RCKWRK, double* WDOT); }

#else

	extern "C" {	void ckwyp_(double* P, double* T, double* Y, int* ICKWRK, double* RCKWRK, double* WDOT); }

#endif

void ChemkinFormationRates(const int ns, const double T, const double P, double* Y, double* R);
