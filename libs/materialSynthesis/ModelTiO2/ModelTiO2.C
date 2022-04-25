/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2022 Alberto Cuoci                                       |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "ModelTiO2.H"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include <iomanip>

const double ModelTiO2::pi_ = std::acos(-1.);		// PI
const double ModelTiO2::kB_ = 1.38064852e-23;		// Boltzmann's constant (m2kg/s2/K)
const double ModelTiO2::Nav_mol_ = 6.02214086e23;	// Avogadro's constant (#/mol)
const double ModelTiO2::Nav_kmol_ = 6.02214086e26;	// Avogadro's constant (#/kmol)

const double ModelTiO2::W_TiO2_ = 79.897999;		// TiO2 molecular weight (kg/kmol or g/mol)
const double ModelTiO2::rhos_ = 4230.;			// solid density (kg/m3) 
const double ModelTiO2::m_TiO2_ = 1.323e-25;		// mass of TiO2 molecule (kg)
const double ModelTiO2::d_TiO2_ = 3.912e-10;		// diameter of TiO2 molecule (m)
const double ModelTiO2::m_TiOH4_ = 1.925e-25;		// mass of TiOH4 molecule (kg)
const double ModelTiO2::d_TiOH4_ = 5.128e-10;		// diameter of TiOH4 molecule (m)

ModelTiO2::ModelTiO2()
{
	// User-defined parameters 
	n0_ = 5;			// number of TiO2 molecules in monomer (default: 5)
	na_ = 0.50;			// nucleation temperature exponent (default: 0.50)
	Ta_ = 0.;			// nucleation activation temperature (default 0 K)
	epsilon_ = 2.2;			// nucleation collision enhancement factor (-)

	// Sub-models
	is_nucleation_ = true;		// nucleation on/off
	is_coagulation_	= true;		// coagulation on/off
	is_sintering_ = true;		// sintering on/off

	// Sintering model: Taus = As*(T^ns)*(dp^4)*exp(Ts_/T);
	As_ = 7.44e16;			// sintering frequency factor (default: 7.44e16 1/s/K/m^(-4))
	ns_ = 1.00;			// sintering temperature exponent (default: 1.00)
	Ts_ = 31000.;			// sintering activation temperature (default 31000 K)

	// Numerical parameters
	N_threshold_ = 1.e3;		// N threshold for calculation of coagulation/sintering contribution (#/m3)
	fv_threshold_ = 1.e-10;		// fv threshold for calculation of properties (-)

	// Precalculation of constant variables dependent on user defined parameters
	Precalculations();
}


void ModelTiO2::Precalculations()
{
	v0_ = n0_*W_TiO2_/Nav_kmol_/rhos_;	// monomer volume (m3)
	d0_ = std::pow(6./pi_*v0_, 1./3.);	// monomer diameter (m)
	s0_ = pi_*d0_*d0_;			// monomer surface area (m2)

	alpha_ = epsilon_*std::sqrt(pi_*kB_/m_TiOH4_)*std::pow(2.*d_TiOH4_,2.);		// (m3/s/K^(1/2))
}


void ModelTiO2::SetNucleationCollisionEnhancementFactor(const double epsilon)
{
	epsilon_ = epsilon;
	Precalculations();
}


void ModelTiO2::SetNucleationActivationTemperature(const double Ta)
{
	Ta_ = Ta;
	Precalculations();
}


void ModelTiO2::SetNucleationTemperatureExponent(const double na)
{
	na_ = na;
	Precalculations();
}

void ModelTiO2::SetSinteringFrequencyFactor(const double As)
{
	As_ = As;
	Precalculations();
}

void ModelTiO2::SetSinteringActivationTemperature(const double Ts)
{
	Ts_ = Ts;
	Precalculations();
}


void ModelTiO2::SetSinteringTemperatureExponent(const double ns)
{
	ns_ = ns;
	Precalculations();
}


void ModelTiO2::SetNumberOfParticles(const unsigned int n)
{
	n0_ = n;
	Precalculations();
}

void ModelTiO2::SetNucleation(const bool flag)
{
	is_nucleation_ = flag;
	Precalculations();
}

void ModelTiO2::SetCoagulation(const bool flag)
{
	is_coagulation_ = flag;
	Precalculations();
}

void ModelTiO2::SetSintering(const bool flag)
{
	is_sintering_ = flag;
	Precalculations();
}

void ModelTiO2::SetNThreshold(const double N_threshold)
{
	N_threshold_ = N_threshold;
	Precalculations();
}

void ModelTiO2::SetFvThreshold(const double fv_threshold)
{
	fv_threshold_ = fv_threshold;
	Precalculations();
}


void ModelTiO2::NucleationSourceTerms(	const double T, const double rho, const double cTiOH4,
					double& r, 
					double& OmegaYs, double& OmegaN, double& OmegaS	)
{
	if (is_nucleation_ == true)
	{
		// Nucleation rate
		const double 	A = 1./4.*alpha_*Nav_kmol_*std::pow(T,0.5-na_)/std::exp(-Ta_/T);	// frequency factor (m3/kmol/s/K^na)
		     		r = A*std::pow(T, na_)*std::exp(-Ta_/T)*(cTiOH4*cTiOH4);		// reaction rate (kmol/m3/s)
		const double QTiO2 = 2.*r;							// formation rate of TiO2 in molar units (kmol/m3/s)
		const double OmegaTiO2 = QTiO2*W_TiO2_;						// formation rate of TiO2 in mass units (kg/m3/s)
		const double I = OmegaTiO2/(rhos_*v0_);						// nucleation rate (1/m3/s)

		// Source terms
		OmegaYs = rhos_/rho*I*v0_;	// (1/s)
		OmegaN = I;			// (1/m3/s)
		OmegaS = I*s0_ ;		// (1/m/s)
	}
	else
	{
		// Source terms
		OmegaYs = 0.;	// mass fraction (1/s)
		OmegaN = 0.;	// particle number density (1/m3/s)
		OmegaS = 0.;	// surface concentration, i.e. surface per unit of volume (1/m/s)
	}
}


void ModelTiO2::CoagulationSourceTerms(	const double T, const double rho, const double mu, 
					const double Ys, const double N, const double S,
					double& OmegaYs, double& OmegaN, double& OmegaS	)
{
	if (N > N_threshold_ && is_coagulation_ == true)
	{
		const double ss = S/N;			// total surface (m2)
		const double fv = rho/rhos_*Ys;		// volume fraction (-) also called total volume concentration (V)
		const double v = fv/N;			// total volume (m3)
	
		// Auxiliary variables
		const double dp = 6.*v/ss;				// primary diameter (m)
		const double np = 6./pi_*v/std::pow(dp,3.);		// number of particles (-)
		const double dc = dp*std::pow(np, 1./1.8);		// collision diameter (m)
		const double D = kB_*T/(3.*pi_*mu*dc);			// particle diffusion coefficient (m2/s)
		const double c = std::sqrt(8.*kB_*T/(pi_*rhos_*v));	// particle velocity (m/s)
		const double L = 8./pi_*D/c;				// mean free path (m)
	
		// Coagulation kernel
		const double g = (std::pow(dc+L,3.)-std::pow(dc*dc+L*L,1.5))/(3.*L*dc)-dc;			// length (m)
		const double Beta = 4.*pi_*dc*D / (0.5*dc/(dc+g*std::sqrt(2.))+D*std::sqrt(2.)/(0.5*c*dc));	// (m3/s)

		// Source terms
		OmegaYs = 0.;			// mass fraction (1/s)
		OmegaN = -0.50*Beta*N*N;	// particle number density (1/m3/s)
		OmegaS = 0.;			// surface concentration, i.e. surface per unit of volume (1/m/s)
	}
	else
	{
		// Source terms
		OmegaYs = 0.;	// mass fraction (1/s)
		OmegaN = 0.;	// particle number density (1/m3/s)
		OmegaS = 0.;	// surface concentration, i.e. surface per unit of volume (1/m/s)
	}
}

void ModelTiO2::SinteringSourceTerms(	const double T, const double rho,
					const double Ys, const double N, const double S,
					double& OmegaYs, double& OmegaN, double& OmegaS	)
{
	if (N > N_threshold_ && is_sintering_ == true)
	{
		const double ss = S/N;					// total surface (m2)
		const double fv = rho/rhos_*Ys;				// volume fraction (-) also called total volume concentration (V)
		const double v = fv/N;					// total volume (m3)
		const double dp = 6.*v/ss;				// primary diameter (m)
		const double ssph = std::pow(v/v0_, 2./3.)*s0_;		// surface area of the completely fused particles (m2)

		// Characteristic sintering time for TiO2 
		// [7.44e16 1/s/K, 1., 31000.] (Bruesser, 2011)
		// [9.75e15 1/s/K, 1., 31000.] (Seto, 1995)
		// [9.11e17 1/s/K, 1., 31000.] (Eggersdorf, 2012)
		const double Taus = As_*std::pow(T,ns_)*std::pow(dp,4.)*std::exp(Ts_/T);	// (1/s)
	
		// Source terms
		OmegaYs = 0.;			// mass fraction (1/s)
		OmegaN = 0.;			// particle number density (1/m3/s)
		OmegaS = -1./Taus*(S-N*ssph);	// surface concentration, i.e. surface per unit of volume (1/m/s)
	}
	else
	{
		// Source terms
		OmegaYs = 0.;	// mass fraction (1/s)
		OmegaN = 0.;	// particle number density (1/m3/s)
		OmegaS = 0.;	// surface concentration, i.e. surface per unit of volume (1/m/s)
	}
}

void ModelTiO2::DiffusionCoefficient(	const double T, const double rho, const double mu, const double Sc,
					const double Ys, const double N, const double S,
					double& Gamma )
{
	Gamma = mu/Sc;

	if (N > N_threshold_)
	{
		const double ss = S/N;			// total surface (m2)
		const double fv = rho/rhos_*Ys;		// volume fraction (-) also called total volume concentration (V)
		const double v = fv/N;			// total volume (m3)
	
		const double dp = 6.*v/ss;				// primary diameter (m)
		const double np = 6./pi_*v/std::pow(dp,3.);		// number of particles (-)
		const double dc = dp*std::pow(np, 1./1.8);		// collision diameter (m)
		const double D = kB_*T/(3.*pi_*mu*dc);			// particle diffusion coefficient (m2/s)

		Gamma = std::max(rho*D, mu/Sc);				// particle diffusion coefficient (kg/m3);
	}

}

void ModelTiO2::Properties(	const double T, const double rho,
				const double Ys, const double N, const double S,
				double& fv, double& dp, double& dc, double& da, double& ssph, double& np, double& Taus )
{
	fv = rho/rhos_*Ys;	// volume fraction (-) also called total volume concentration (V)

	if (fv > fv_threshold_)
	{
		const double ss = S/N;						// total surface (m2)   
		const double v = fv/N;						// total volume (m3)
	
		dp = 6.*v/ss;							// primary diameter (m)
		da = std::pow(6.*v/(pi_*N), 1./3.);				// aggregate diameter (m)
		np = 6./pi_*v/std::pow(dp,3.);					// number of particles (-)
		dc = dp*std::pow(np, 1./1.8);					// collision diameter (m)
		ssph = std::pow(v/v0_, 2./3.)*s0_;				// surface area of the completely fused particles (m2)
		Taus = As_*std::pow(T,ns_)*std::pow(dp,4.)*std::exp(Ts_/T);	// (1/s)
	}
	else
	{
		dp = 0.;	// primary diameter (m)
		da =  0.;	// aggregate diameter (m)
		np = 0.;	// number of particles (-)
		dc = 0.;	// collision diameter (m)
		ssph = 0.;	// surface area of the completely fused particles (m2)
		Taus = 0.;	// (1/s)
	}
}

void ModelTiO2::SummaryOnScreen()
{	
	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << "                               TiO2 Model Summary                                         " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << " * TiO2 molecular weight (kg/kmol): " << W_TiO2_ << std::endl;
	std::cout << " * TiO2 mass (kg): " << m_TiO2_ << std::endl;
	std::cout << " * TiO2 diameter (m): " << d_TiO2_ << std::endl;
	std::cout << " * TiOH4 mass (kg): " << m_TiOH4_ << std::endl;
	std::cout << " * TiOH4 diameter (m): " << d_TiOH4_ << std::endl;
	std::cout << " * Solid density (kg/m3): " << rhos_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Processes" << std::endl;
	std::cout << "    + Nucleation: " << is_nucleation_ << std::endl;
	std::cout << "    + Coagulation: " << is_coagulation_ << std::endl;
	std::cout << "    + Sintering: " << is_sintering_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Monomer data" << std::endl;
	std::cout << "    + Number of molecules: " << n0_ << std::endl;
	std::cout << "    + Diameter (m): " << d0_ << std::endl;
	std::cout << "    + Surface (m2): " << s0_ << std::endl;
	std::cout << "    + Volume (m3): " << v0_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Nucleation parameters" << std::endl;
	std::cout << "    + Collision enhancement factor (-): " << epsilon_ << std::endl;
	std::cout << "    + Temperature exponent (-): " << na_ << std::endl;
	std::cout << "    + Activation temperature (K): " << Ta_ << std::endl;
	std::cout << "    + Nucleation coefficient (m3/s/K^(1/2)): " << alpha_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Coagulation parameters" << std::endl;
	std::cout << "    + Frequency factor (s,K): " << As_ << std::endl;
	std::cout << "    + Temperature exponent: " << ns_ << std::endl;
	std::cout << "    + Activation temperature (K): " << Ts_ << std::endl;
	std::cout << std::endl;
	std::cout << " * Numerical parameters" << std::endl;
	std::cout << "    + Threshold on N (#/m3): " << N_threshold_ << std::endl;
	std::cout << "    + Threshold on fv (-): " << fv_threshold_ << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl; 
	std::cout << std::endl;
}

void ModelTiO2::FatalErrorMessage(const std::string message)
{
	std::cout << "ModelTiO2: Fatal error" << std::endl;
	std::cout << "Error message: " << message << std::endl;
	std::cout << "Press enter to exit..." << std::endl;
	getchar();
	exit(-1);
}
